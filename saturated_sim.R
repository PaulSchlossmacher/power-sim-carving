# Clear all variables
rm(list = ls())


#Local, user specific path, that should work for both of us:
Local_path<-getwd()
hdi_adjustments_path<-paste(Local_path, "/Multicarving-Christoph/inference/hdi_adjustments.R", sep="")
carving_path<-paste(Local_path, "/Multicarving-Christoph/inference/carving.R", sep="")
sample_from_truncated_path<-paste(Local_path, "/Multicarving-Christoph/inference/sample_from_truncated.R", sep="")
tryCatchWE_path<-paste(Local_path, "/Multicarving-Christoph/inference/tryCatch-W-E.R", sep="")

#Different paths here, because they're "our own" functions
SNTN_distribution_path<-paste(Local_path, "/SNTN_distribution.R", sep="")
split_select_function_path<-paste(Local_path, "/split_select.R", sep="")
carve_linear_path<-paste(Local_path, "/carve_linear.R", sep="")


library(MASS)
library(mvtnorm)
library(glmnet)
library(Matrix)
library(tictoc)
library(hdi)
library(selectiveInference)
library(doSNOW)
library(parallel)
library(doParallel)
library(doRNG)
library(truncnorm)
library(git2r)
library(ggplot2)
library(Rmpfr)

source(hdi_adjustments_path)
source(carving_path)
source(sample_from_truncated_path)
source(tryCatchWE_path)
source(SNTN_distribution_path)
source(split_select_function_path)
source(carve_linear_path)




carve.saturated <- function(x, y, B = 50, fraction = 0.9, gamma = ((1:B)/B)[((1:B)/B) >= 0.05], FWER = FALSE, ci.level = 0.95,
                                     family = "gaussian", model.selector = lasso.cvcoef,
                                     args.model.selector = list(intercept = TRUE, standardize = FALSE),
                                     se.estimator = "modwise", args.se.estimator = list(df.corr = TRUE, intercept = TRUE, standardize = FALSE),
                                     args.lasso.inference = list(sigma = NA), ci.timeout = 10, split.pval = TRUE,
                                     classical.fit = lm.pval, args.classical.fit = NULL, args.classical.ci = NULL,
                                     parallel = FALSE, ncores = getOption("mc.cores", 2L),
                                     return.nonaggr = FALSE, return.selmodels = FALSE, verbose = FALSE) {
  # routine to split the data, select a model, calculate carving p-values B times and determine the corresponding CI
  # Input
  # x (matrix): matrix of predictors
  # y (vector): response vector
  # B (integer): number of splits
  # fraction (numeric in (0,1)): fraction used for selection
  # gamma (numeric or vector of numeric in (0,1]): quantiles to consider, if several, additional penalty is applied
  # FWER (boolean): shall a FWER correction be applied
  # ci.level (numeric in (0,1)): level of the confindence interval
  # family (string): "gaussian" or "binomial"
  # model.selector (function): how the model is chosen (must be some version of Lasso)
  # args.model.selector (list): additional arguments for the selection process
  # se.estimator (string): how sigma is estimated, "1se", "modwise", "min" or "None"
  # args.se.estimator (list): additional arguments to estimate sigma
  # args.lasso.inference (list): additional arguments for inference after Lasso
  # ci.timeout (numeric > 0): maximum time to search for an uncovered point before setting the bound to Inf/ -Inf
  # split.pval (boolean): shall p-values and confidence intervals for splitting be determined as well
  # classical.fit (function): function to calculate splitting p-values
  # args.classical.fit (list): additional arguments for calculating splitting p-values
  # args.classical.ci (list): additional arguments for calculating splitting CI
  # parallel (boolean): whether to parallelize the splits (CI calculation is never parallelized)
  # ncores (integer): number of cores for parallelization
  # return.nonaggr (boolean): shall raw p-values be returned
  # return.sel.models (boolean): shall the information, which model was selected be returned
  # verbose (boolean): whether to print key steps
  # Output (if split.pval = TRUE,  a list of two output elements is created)
  # pval.corr (p - vector): multicarving / multisplitting p-values
  # gamma.min (p - vector): which value of gamma minimized the p-value
  # ci.level (numeric in (0,1)): level of the confidence interval
  # lci (p - vector): lower end of confidence interval
  # uci (p - vector): upper end of confidence interval
  # optional: pvals.nonaggr (B x p matrix): p-values before aggregation
  # optional: sel.models (boolean B x p matrix): TRUE if variable was selected in given split
  # FWER (boolean): was a multiplicity correction applied?
  # only for carving: vlo (B x p matrix): lower end of constrained region
  # only for carving: vup (B x p matrix): upper end of constrained region
  # only for carving: centers (B x p matrix): estimate of parameter in given model
  # ses (B x p matrix): estimate of standard error of parameter in given model
  # method ("multi.carve" or "multi.split"): additional output information
  # call: function call to obtain this carving result
  

  #---------- Note: I'm simplifying this code under the assumption that split.pval=FALSE, sigma is provided
  
  args.model.selector$family <- family
  args.lasso.inference$family <- family
  
  # provided sigma has priority over se estimator
  se.estimator <- "None"
  globalSigma <- args.lasso.inference$sigma

  
  n <- nrow(x)
  p <- ncol(x)
  n.left <- floor(n * fraction)
  n.right <- n - n.left
  stopifnot(n.left >= 1, n.right >= 0)
  oneSplit <- function(b) {
    
    pvals.v <- rep(1, p)
    
    sel.models <- logical(p)
    vlo.v <- rep(-Inf, p)
    vup.v <- rep(Inf, p)
    estimates.v <- rep(NA, p)
    sescarve.v <- rep(Inf, p)
    lci.v <- rep(-Inf, p)
    uci.v <- rep(Inf, p)
    centers.v <- rep(NA, p)
    ses.v <- rep(Inf, p)
    df.res <- NA
    try.again <- TRUE
    thresh.count <- 0L
    threshn <- 1e-7
    continue <- TRUE
    split.again <- TRUE
    split.count <- 0
    while (split.again) {
      split.again <- FALSE
      split <- sample.int(n, size = n.left)
      x.left <- x[split, ]
      y.left <- y[split]
      x.right <- x[-split, ]
      y.right <- y[-split]
      
      output <- do.call(model.selector, args = c(list(x = x.left, 
                                                      y = y.left), args.model.selector))
      sel.model <- output$sel.model
      beta <- output$beta
      lambda <- output$lambda
      
      fit.again <- TRUE
      thresh.count <- 0
      p.sel <- length(sel.model)
      if (p.sel == 0) fit.again <- FALSE
      
      while (fit.again) {
        fit.again <- FALSE
        checktry <- tryCatch_W_E(constraint.checker(x.left, y.left, beta, 0, lambda, family,
                                                    intercept = args.model.selector$intercept), TRUE)
        if (is.null(checktry$error)) {
          check <- checktry$value
        } else {
          check <- TRUE
          split.again <- TRUE
          warning(paste(checktry$error, p.sel, " variables selected with ",
                        n.left, "data points, splitting again"))
        }
        if (!check) {
          if (verbose) 
            cat("......fit again...\n")
          fit.again <- TRUE
          thresh.count <- thresh.count + 1
          if (thresh.count > 2) {
            warning("Giving up reducing threshhold")
            break()
          }
          threshn <- 1e-7 / (100) ^ thresh.count
          fit <- glmnet(x = x.left,y = y.left,standardize = args.model.selector$standardize,
                        intercept = args.model.selector$intercept, thresh = threshn, family = family)
          if (verbose) cat(threshn,"\n")
          coefs <- coef(fit, x = x.left, y = y.left, s = lambda / n.left, exact = TRUE,
                        standardize = args.model.selector$standardize,
                        intercept = args.model.selector$intercept, thresh = threshn, family = family)
          beta <- coefs[-1]
          sel.model <- which(abs(beta) > 0)
          
          p.sel <- length(sel.model)
          if (p.sel == 0) fit.again <- FALSE
          warning(paste("reducing threshold", thresh.count, "to", threshn, sep = " "))
        }
      }
      
      p.sel <- length(sel.model)
      # use new split in case of singularity. This is mostly an issue for discrete x.
      if (args.model.selector$intercept) {
        if ((p.sel > 0 && (rankMatrix(cbind(rep(1, n.left), x.left[, sel.model]))[[1]]< (p.sel + 1) ||
                           (p.sel < n.right - 1 && rankMatrix(cbind(rep(1, n.right), x.right[, sel.model]))[[1]] < (p.sel + 1)))) ||
            fit.again) split.again <- TRUE
      } else {
        if ((p.sel > 1 && (rankMatrix(x.left[, sel.model])[[1]] < (p.sel) ||
                           (p.sel < n.right  && rankMatrix(x.right[, sel.model])[[1]]< (p.sel)))) ||
            fit.again) split.again <- TRUE
      }
      if (split.again) {
        reason <- character(0)
        if (args.model.selector$intercept){
          if (rankMatrix(cbind(rep(1, n.left), x.left[,sel.model]))[[1]] < (p.sel + 1)) reason <- c(reason, "x1 rank")
          if (p.sel < n.right - 1 && rankMatrix(cbind(rep(1, n.right), x.right[,sel.model]))[[1]] < (p.sel + 1)) reason <- c(reason, "x2 rank")
        } else {
          if (rankMatrix( x.left[,sel.model])[[1]] < (p.sel)) reason <- c(reason, "x1 rank")
          if (p.sel < n.right && rankMatrix(x.right[,sel.model])[[1]] < (p.sel)) reason <- c(reason, "x2 rank")
        }
        
        if (fit.again) reason <- c(reason, "fit")
        if (!is.null(checktry$error)) reason <- c(reason, "error while checking")
        split.count <- split.count + 1
        if (split.count > 4) {
          stop(paste("More than 5 splits needed, final reason:", reason))
        }
        if (verbose) 
          cat("......Splitting again...\n")
        warning(paste("Splitting again ", split.count, "reason", reason))
      }
    }
    
    if (p.sel > 0) {
      fLItry <- tryCatch_W_E(do.call(carve.lasso,
                                     args = c(list(X = x, y = y,ind = split, beta = beta, tol.beta = 0,
                                                   lambda = lambda, intercept = args.model.selector$intercept,
                                                   selected = FALSE), args.lasso.inference)), 0)
      if (!is.null(fLItry$error)) {
        warning(paste("Failed to infer a split, due to:", fLItry$error, sep=" "))
        pvals.v[] = NA
        list(pvals = pvals.v, sel.models = sel.models, centers = centers.v, 
             ses = ses.v, df.res = df.res, lci = lci.v, uci = uci.v, sescarve = sescarve.v,
             vlo = vlo.v, vup = vup.v, estimates= estimates.v, split = split)
      } else if (!is.null(fLItry$warning)) {
        for (war in unique(fLItry$warning)) {
          warning(paste("Split", b, ":", war, sep = " "))
        }
      }
      fLI <- fLItry$value
      sel.pval1 <- fLI$pv
      sel.vlo <- fLI$vlo
      sel.vup <- fLI$vup
      sel.sescarve <- fLI$ses
      sel.estimates <- fLI$estimates
      
      if (any(is.na(sel.pval1))) {
        stop("The carve procedure returned a p-value NA")
      } 
      
      if (length(sel.pval1) != p.sel) { 
        stop(paste("The carve procedure didn't return the correct number of p-values for the provided submodel.",
                   p.sel, length(sel.pval1)))
      }
      if (!all(sel.pval1 >= 0 & sel.pval1 <= 1)) {
        stop("The carve procedure returned values below 0 or above 1 as p-values")
      }
      sel.pval1 <- 2 * pmin(sel.pval1, 1 - sel.pval1)
      if (FWER) {
        sel.pval1 <- pmin(sel.pval1 * p.sel, 1) # for FWER
      } else {
        sel.pval1 <- pmin(sel.pval1, 1) # for FCR
      }
      
      pvals.v[sel.model] <- sel.pval1
      vlo.v[sel.model] <- sel.vlo
      vup.v[sel.model] <- sel.vup
      estimates.v[sel.model] <- sel.estimates
      sescarve.v[sel.model] <- sel.sescarve
      
      if (return.selmodels) 
        sel.models[sel.model] <- TRUE
    }
    
    if (p.sel == 0) {
      if (verbose) 
        cat("......Empty model selected. That's ok...\n")
    }
    list(pvals = pvals.v, sel.models = sel.models, centers = centers.v, 
         ses = ses.v, df.res = df.res, lci = lci.v, uci = uci.v, sescarve = sescarve.v,
         vlo = vlo.v, vup = vup.v, estimates= estimates.v, split = split)
  }
  split.out <- if (parallel) {
    stopifnot(isTRUE(is.finite(ncores)), ncores >= 1L)
    if (verbose) 
      cat("...starting parallelization of sample-splits for selection and constraints\n")
    mclapply(1:B, oneSplit, mc.cores = ncores)
  } else {
    if (verbose)
      cat("...selecting models and determing constraints\n")
    lapply(1:B, oneSplit)
  }
  myExtract <- function(name) {
    matrix(unlist(lapply(split.out, "[[", name)), nrow = B, 
           byrow = TRUE)
  }
  if (verbose) 
    cat("...determing confidence intervals\n")

    pvals <- myExtract("pvals")
  colnames(pvals) <- colnames(x)
  if (return.selmodels) {
    sel.models <- myExtract("sel.models")
    colnames(sel.models) <- colnames(x)
  } else {
    sel.models <- NA
  }
  vlo <- myExtract("vlo")
  vup <- myExtract("vup")
  estimates <- myExtract("estimates")
  sescarve <- myExtract("sescarve")
  pvals.current <- which.gamma <- numeric(p)
  for (j in 1:p) {
    quant.gamma <- quantile(pvals[, j], gamma, na.rm = TRUE,type = 3) / gamma
    penalty <- if (length(gamma) > 1) 
      (1 - log(min(gamma)))
    else 1
    pvals.pre <- min(quant.gamma) * penalty
    pvals.current[j] <- pmin(pvals.pre, 1)
    which.gamma[j] <- which.min(quant.gamma)
  }
  names(pvals.current) <- names(which.gamma) <- colnames(x)
  vars <- ncol(vlo)
  s0 <- if (any(is.na(sel.models))) 
    NA
  else apply(sel.models, 1, sum)
  new.ci <- mapply(aggregate.ci.saturated, vlo = split(vlo, rep(1:vars, each = B)),
                   vup = split(vup, rep(1:vars, each = B)), 
                   centers = split(estimates, rep(1:vars, each = B)), 
                   ses = split(sescarve, rep(1:ncol(sescarve), each = B)), 
                   gamma.min = min(gamma), multi.corr = FALSE, verbose = FALSE, timeout = ci.timeout,
                   s0 = list(s0 = s0), ci.level = ci.level, var = 1:vars)
  lci.current <- t(new.ci)[, 1]
  uci.current <- t(new.ci)[, 2]
  if (!return.nonaggr) 
    pvals <- NA
  names(lci.current) <- names(uci.current) <- names(pvals.current)
  if (return.selmodels) {
    keep <- c("return.selmodels", "x", "y", "gamma", "split.out", "FWER",
              "pvals", "pvals.current", "which.gamma", "sel.models", "ci.level", 
              "lci.current", "uci.current", "vlo", "vup", "sescarve", "estimates")
    rm(list = setdiff(names(environment()), keep))
  }
  # structure(list(pval.corr = pvals.current, gamma.min = gamma[which.gamma], 
  #                ci.level = ci.level, lci = lci.current, uci = uci.current,
  #                pvals.nonaggr = pvals, sel.models = sel.models, FWER = FWER,
  #                vlo = vlo, vup = vup, centers = estimates, ses = sescarve,
  #                method = "multi.carve", call = match.call()), class = "carve")
  # 
  return(list(pval.corr = pvals.current, gamma.min = gamma[which.gamma], 
                 ci.level = ci.level, lci = lci.current, uci = uci.current,
                 pvals.nonaggr = pvals, sel.models = sel.models, FWER = FWER,
                 vlo = vlo, vup = vup, centers = estimates, ses = sescarve,
                 method = "multi.carve", call = match.call()), class = "carve")
  }
}

# toeplitz
n <- 100
p <- 200
rho <- 0.6
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
sel.index <- c(1, 5, 10, 15, 20)
ind <- sel.index
beta <- rep(0, p)
beta[sel.index] <- 1
sparsity <- 5
set.seed(42) # to make different methods comparable, fix the x-matrix
x <- mvrnorm(n, rep(0, p), Cov) 
print (x[1,1])
# should create the right x on D-MATH server, x[1 ,1] = 0.958
sigma <- 2
y.true <- x %*% beta
y <- y.true + sigma * rnorm(n)


B=1
frac=0.7

#gammavec <- round(seq(ceiling(0.05 * B) / B, 1 - 1 / B, by = 1 / B), 2)
#+gammavec[1] <- 0.05999 # due to inconsistency for multisplitting CI

gammavec <- 0.05999

results<-multi.carve.ci.saturated(x, y, B = B, fraction = frac, ci.level = 0.95, 
                         model.selector = lasso.cvcoef, classical.fit = lm.pval,
                         parallel = FALSE, ncores = getOption("mc.cores", 2L), gamma = gammavec,
                         args.model.selector = list(standardize = FALSE, intercept = TRUE,
                                                    tol.beta = 0, use.lambda.min = FALSE),
                         verbose = FALSE, ci.timeout = 10, FWER = FALSE, split.pval = TRUE,
                         return.selmodels = TRUE, return.nonaggr = TRUE)

results
