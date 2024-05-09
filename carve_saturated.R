# carve.saturated

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
  
  if (!(se.estimator %in% c("1se", "min", "modwise", "None")))
    stop("Sigma estimator must be one of \"1se\", \"min\", \"modwise\" or \"None\" ")
  
  if (!(family %in% c("gaussian", "binomial")))
    stop ("Invalid family provided, can only deal with gaussian and binomial")
  
  args.model.selector$family <- family
  args.lasso.inference$family <- family
  
  if (is.null(args.lasso.inference$sigma)) args.lasso.inference$sigma <- NA
  if (family == "gaussian"){
    if (se.estimator == "None" && is.na(args.lasso.inference$sigma)) stop("Neither SE estimator type nor sigma provided for Gaussian family. This is not ok")
    if (is.na(args.lasso.inference$sigma)) {
      if (se.estimator %in% c("1se", "modwise")){
        use.lambda.min = FALSE
      } else {
        use.lambda.min = TRUE
      }
      estSigma <- do.call(estimateSigma.flex,
                          args = c(list(x = x, y = y, use.lambda.min = use.lambda.min), args.se.estimator))
      globalSigma <- estSigma$sigmahat
      args.lasso.inference$sigma <- globalSigma
    } else {
      # provided sigma has priority over se estimator
      se.estimator <- "None"
      globalSigma <- args.lasso.inference$sigma
    }
  }
  
  n <- nrow(x)
  p <- ncol(x)
  n.left <- floor(n * fraction)
  n.right <- n - n.left
  stopifnot(n.left >= 1, n.right >= 0)
  oneSplit <- function(b) {
    if (verbose) 
      cat("...split", b, "\n")
    if (split.pval) {
      pvals.v <- matrix(1, nrow = 2, ncol = p)
    } else {
      pvals.v <- rep(1, p)
    }
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
          
          if (family == "binomial") beta <- coefs
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
    
    if (se.estimator == "modwise" && family == "gaussian") {
      if (length(beta) == p + 1) beta <- beta[-1]
      if (args.model.selector$intercept){
        RSS <- sum((scale(y, T, F) - scale(x, T, F) %*% beta) ^ 2)
        if (args.se.estimator$df.corr) {
          den <- n - p.sel - 1
        } else {
          den <- n
        }
        sigma.model <- sqrt(RSS / den)
      } else {
        RSS <- sum((y- x %*% beta) ^ 2)
        if (args.se.estimator$df.corr) {
          den <- n - p.sel
        } else {
          den <- n
        }
        sigma.model <- sqrt(RSS / den)
      }
      estSigma <- sigma.model
      args.lasso.inference$sigma <- sigma.model
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
      if (split.pval) {
        x.right <- x[-split, ]
        y.right <- y[-split]
        if (args.model.selector$intercept) {
          bound <- n.right-1
        } else {
          bound <- n.right
        }
        if (p.sel < bound) {
          sel.pval2try <- tryCatch_W_E(do.call(classical.fit, 
                                               args = c(list(x = x.right[, sel.model], y = y.right), args.classical.fit)),
                                       rep(NA, p.sel))
          sel.pval2 <- sel.pval2try$value
          if (!is.null(sel.pval2try$error)) {
            warning(paste(sel.pval2try$error, "while caluclatng split p-values", sep=" "))
          }
          NAs <- FALSE
          if (any(is.na(sel.pval2))) NAs <- TRUE
          # do not stop if splitting leads to NA
          if (length(sel.pval2) != p.sel) 
            stop("The classical.fit function didn't return the correct number of p-values for the provided submodel.")
          if (!all(sel.pval2 >= 0 & sel.pval2 <= 1) && !NAs) 
            stop("The classical.fit function returned values below 0 or above 1 as p-values")
          if (FWER) {
            sel.pval2 <- pmin(sel.pval2 * p.sel, 1) # for FWER
          } else {
            sel.pval2 <- pmin(sel.pval2, 1) # for FCR
          }
          tmp.fit.lm <- lm(y.right ~ x.right[, sel.model], 
                           args.classical.fit)
          a <- (1 - ci.level)/2
          a <- c(a, 1 - a)
          fac <- qt(a, tmp.fit.lm$df.residual)
          sel.ses <- sqrt(diag(vcov(tmp.fit.lm)))[-1]
          sel.centers <- coef(tmp.fit.lm)[-1]
          sel.ci <- sel.centers + sel.ses %o% fac
          centers.v[sel.model] <- sel.centers
          lci.v[sel.model] <- sel.ci[, 1]
          uci.v[sel.model] <- sel.ci[, 2]
          ses.v[sel.model] <- sel.ses
          df.res <- tmp.fit.lm$df.residual
          pvals.v[1, sel.model] <- sel.pval1
          pvals.v[2, sel.model] <- sel.pval2
          vlo.v[sel.model] <- sel.vlo
          vup.v[sel.model] <- sel.vup
          estimates.v[sel.model] <- sel.estimates
          sescarve.v[sel.model] <- sel.sescarve
        } else {
          pvals.v[1, sel.model] <- sel.pval1
          vlo.v[sel.model] <- sel.vlo
          vup.v[sel.model] <- sel.vup
          estimates.v[sel.model] <- sel.estimates
          sescarve.v[sel.model] <- sel.sescarve
        }
      } else {
        pvals.v[sel.model] <- sel.pval1
        vlo.v[sel.model] <- sel.vlo
        vup.v[sel.model] <- sel.vup
        estimates.v[sel.model] <- sel.estimates
        sescarve.v[sel.model] <- sel.sescarve
      }
      
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
  if (split.pval) { 
    ls <- list()
    pvalsall <- array(unlist(lapply(split.out, "[[", "pvals")), dim = c(2, p, B))
    for (icf in  1:2) {
      pvals <- t(pvalsall[icf, , ])
      colnames(pvals) <- colnames(x)
      if (return.selmodels) {
        sel.models <- myExtract("sel.models")
        colnames(sel.models) <- colnames(x)
      } else {
        sel.models <- NA
      }
      if (icf == 1) {
        vlo <- myExtract("vlo")
        vars <- ncol(vlo)
        vup <- myExtract("vup")
        sescarve <- myExtract("sescarve")
        estimates <- myExtract("estimates")
      } else {
        lci <- myExtract("lci")
        uci <- myExtract("uci")
        centers <- myExtract("centers")
        ses <- myExtract("ses")
      }
      
      df.res <- unlist(lapply(split.out, `[[`, "df.res"))
      pvals.current <- which.gamma <- numeric(p)
      for (j in 1:p) {
        if (any(!is.na(pvals[, j]))) {
          quant.gamma <- quantile(pvals[, j], gamma, na.rm = TRUE, type = 3) / gamma
          penalty <- if (length(gamma) > 1)
            (1 - log(min(gamma)))
          else 1
          pvals.pre <- min(quant.gamma) * penalty
          pvals.current[j] <- pmin(pvals.pre, 1)
          which.gamma[j] <- which.min(quant.gamma)
        } else {
          pvals.current[j] <- NA
          which.gamma[j] <- NA
        }
      }
      names(pvals.current) <- names(which.gamma) <- colnames(x)
      s0 <- if (any(is.na(sel.models)))
        NA
      else apply(sel.models, 1, sum)
      if (icf == 1) {
        new.ci <- mapply(aggregate.ci.saturated, vlo = split(vlo, rep(1:vars, each = B)),
                         vup = split(vup, rep(1:vars, each = B)), 
                         centers = split(estimates, rep(1:vars, each = B)), 
                         ses = split(sescarve, rep(1:ncol(sescarve), each = B)), 
                         gamma.min = min(gamma), multi.corr = FWER, verbose = FALSE, timeout = ci.timeout,
                         s0 = list(s0 = s0), ci.level = ci.level, var = 1:vars)
      } else {
        new.ci <- mapply(hdi:::aggregate.ci, lci = split(lci, rep(1:vars, each = B)),
                         rci = split(uci, rep(1:vars, each = B)),
                         centers = split(centers, rep(1:vars, each = B)),
                         ses = split(ses, rep(1:ncol(ses), each = B)), df.res = list(df.res = df.res),
                         gamma.min = min(gamma), multi.corr = FALSE, verbose = FALSE,
                         s0 = list(s0 = s0), ci.level = ci.level, var = 1:vars)
      }
      lci.current <- t(new.ci)[, 1]
      uci.current <- t(new.ci)[, 2]
      names(lci.current) <- names(uci.current) <- names(pvals.current)
      
      if (!return.nonaggr)
        pvals <- NA
      if (return.selmodels) {
        if (icf == 2) {
          keep <- c("return.selmodels", "x", "y", "gamma", "split.out",
                    "pvals", "pvals.current", "which.gamma", "sel.models",
                    "ls", "icf", "ci.level", "estimates", "FWER",
                    "lci.current", "uci.current", "vlo", "vup", "sescarve", "ses")
          rm(list = setdiff(names(environment()), keep))
        }
      }
      if (icf == 1) {
        ls[[icf]] <- structure(list(pval.corr = pvals.current, gamma.min = gamma[which.gamma], 
                                    ci.level = ci.level, lci = lci.current, uci = uci.current,
                                    pvals.nonaggr = pvals, sel.models = sel.models, FWER = FWER,
                                    vlo = vlo, vup = vup, centers = estimates, ses = sescarve,
                                    method = "multi.carve", call = match.call()), class = "carve")
      } else {
        ls[[icf]] <- structure(list(pval.corr = pvals.current, gamma.min = gamma[which.gamma],
                                    ci.level = ci.level, lci =  lci.current, uci = uci.current,
                                    pvals.nonaggr = pvals, sel.models = sel.models, FWER = FWER,
                                    ses = ses, method = "multi.split", call = match.call()), class = "carve")
      }
    }
    return(ls)
  } else {
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
    structure(list(pval.corr = pvals.current, gamma.min = gamma[which.gamma], 
                   ci.level = ci.level, lci = lci.current, uci = uci.current,
                   pvals.nonaggr = pvals, sel.models = sel.models, FWER = FWER,
                   vlo = vlo, vup = vup, centers = estimates, ses = sescarve,
                   method = "multi.carve", call = match.call()), class = "carve")
  }
}