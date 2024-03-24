
#This is actually the one_split function from Christophs multi.carve function in hdi_adjustments

split.select <- function(x,y,fraction = 0.9, family = "gaussian",model.selector = lasso.cvcoef, 
                            args.model.selector = list(intercept = FALSE, standardize = FALSE)) {
  n <- nrow(x)
  p <- ncol(x)
  n.left <- floor(n * fraction)
  n.right <- n - n.left
  stopifnot(n.left >= 1, n.right >= 0)
  sel.models <- logical(p)
  try.again <- TRUE
  split.count <- 0
  thresh.count <- 0L
  threshn <- 1e-7
  continue <- TRUE
  split.again <- TRUE
  while (split.again) {
    split.again <- FALSE
    #choose random indices to split observations
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
    # for empty model, active constraints are trivially fulfilled
    if (p.sel == 0) fit.again <- FALSE
    
    while (fit.again) {
      fit.again <- FALSE
      checktry <- tryCatch_W_E(constraint.checker(x.left, y.left, beta, 0, lambda,
                                                  family, intercept = args.model.selector$intercept), TRUE)
      if (is.null(checktry$error)) {
        check <- checktry$value
      } else {
        # if run into numerical instability when checking constraints
        check <- TRUE
        split.again <- TRUE
        warning(paste(checktry$error, p.sel, " variables selected with ",
                      n.left, "data points, splitting again"))
      }
      if (!check) {
        cat("......fit again...\n")
        fit.again <- TRUE
        thresh.count <- thresh.count + 1
        if (thresh.count > 2) {
          warning("Giving up reducing threshhold")
          break()
        }
        threshn <- 1e-7 / (100) ^ thresh.count
        fit <- glmnet(x = x.left, y = y.left, standardize = args.model.selector$standardize,
                      intercept = args.model.selector$intercept, thresh = threshn,family = family)
        cat(threshn,"\n")
        coefs <- coef(fit,x = x.left,y = y.left,s = lambda/n.left,exact = TRUE,
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
      cat("......Splitting again...\n")
      warning(paste("Splitting again ", split.count, "reason", reason))
    }
  }
  sel.models[sel.model] <- TRUE
  return (list(sel.models = sel.models, split = split, beta = beta, lambda = lambda))
}

