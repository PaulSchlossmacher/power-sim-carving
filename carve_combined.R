
#' splits the data, calculates p-values of carving estimator as in Drydale's paper
#'
#' @param x (matrix dim: n x p) the full design matrix
#' @param y (vector dim: n) the full response vector
#' @param split (vector dim: n_a) indices of observations used for selection
#' @param beta (vector dim: p) full beta obtained from lasso selection
#' @param lambda (numeric) lambda from lasso selection
#' @param sigma_squ (numeric) true variance of data generating process y~N(0,sigma^2*I_n)
#'
#' @return list containing p-values of combined carving estimator and all of the inputs used for SNTN_distribution

carve.comb <- function(x, y, split, beta, lambda,sigma_squ){
  
  #Split the data
  n <- length(y)
  p <- length(beta)
  n.a <- length(split)
  n.b <- n-n.a
  x.a <- x[split, ]
  y.a <- y[split, ]
  x.b <- matrix(x[-split, ], nrow = n.b)#matrix for the edge case when n.b = 1
  y.b <- y[-split, ]

  #Sigma gets chosen in accordance with the distribution of y_A in Lemma 2.4.1
  Sigma <- diag(n.a)*sigma_squ
  c1<-n.b/n
  c2<-n.a/n
  
  chosen <-  which(abs(beta) > 0) # selected variables
  s <- length(chosen)
  b.signs <- sign(beta[chosen])
  
  if (s == 0)
    stop("0 variables were chosen by the Lasso")
  
  #extract active variables from both splits
  x.Ma <- x.a[, chosen]#(n.a x s)
  x.Mb <- x.b[, chosen]#(n.b x s)
  
  #extract inactive variables from both splits
  x_Ma <- x.a[, -chosen]#(n.a x (p-s))
  x_Mb <- x.b[, -chosen]
  
  #compute the moore penrose inverse of active variables in both splits
  x.Ma.i <- ginv(x.Ma)#(s x na)
  x.Ma.ti <- ginv(t(x.Ma))
  x.Mb.i <- ginv(x.Mb)
  
  #compute the projection matrix onto active columns
  p.Ma <- x.Ma %*% x.Ma.i#(n.a x n.a)
  
  #compute the beta_carve from drysdales paper
  beta_split<-x.Mb.i %*% y.b
  beta_posi<-x.Ma.i %*% y.a
  
  beta_comb <- ((n.a/n)*x.Ma.i %*% y.a + (n.b/n)*x.Mb.i %*% y.b)
  
  #Get inactive affine constraints on split A, as this is the group where we are doing POSI(Lee et al. page 8)
  #This is commented as it is not used further. For more details check the paper
  # A.0.up <- (t(x_Ma) %*% (diag(n.a) - p.Ma))#(p-s) x n.a
  # A.0 <- 1/lambda*rbind(A.0.up, -A.0.up) #2*(p-s) x n.a
  # b.0.temp <-t(x_Ma) %*% x.Ma.ti %*% b.signs
  # b.0.up <- rep(1,p-s) - b.0.temp
  # b.0.lo <- rep(1,p-s) + b.0.temp
  # b.0 <- rbind(b.0.up, b.0.lo)#2*(p-s) x n.a
  
  #Get active affine constraints on split A
  C <- solve(crossprod(x.Ma,x.Ma))#(X_Ma^T*X_Ma)^-1
  A.1 <- -diag(x = b.signs, nrow = s) %*% x.Ma.i # (s x n.a)
  b.1 <- -lambda*diag(x = b.signs, nrow = s) %*% C %*% b.signs
  
  #A <- rbind(A.0, A.1)# (2p-s) x n.a
  #b <- rbind(b.0,b.1)
  A <- A.1
  b <- b.1
  
  #Calculating truncation limits
  vup <- rep(0,s)
  vlo <- rep(0,s)
  for (i in 1:s){
   
    eta <- x.Ma.i[i,]
    c <- (Sigma %*% eta) / as.numeric((t(eta) %*% Sigma %*% eta))
    z <- (diag(n.a) - c %*% t(eta)) %*% y.a
    den <- A%*%c
    resid <- b-A %*% z
    #We do not consider the set V^0(z) as defined in Lee p.10
    ind.vup <- (den > 0)
    ind.vlo <- (den < 0)
    
    if (any(ind.vup)){
      vup[i] <- min(resid[ind.vup]/den[ind.vup])
    }else {
      vup[i] <- Inf
    }
    
    if (any(ind.vlo)){
      vlo[i] <- max(resid[ind.vlo]/den[ind.vlo])
    }else {
      vlo[i] <- -Inf
    }
  }
  
  tau.M <- sigma_squ
  
  tau.1 <- diag(tau.M*solve(t(x.Mb)%*%x.Mb))
  tau.2 <- diag(tau.M*solve(t(x.Ma)%*%x.Ma))
  
  #Null hypothesis
  theta.1 <- rep(0,s)
  theta.2 <- rep(0,s)
  
  cdf<-SNTN_CDF(z=beta_comb,
                mu1=theta.1,
                tau1=tau.1,
                mu2=theta.2,
                tau2=tau.2,
                a=vlo,
                b=vup,
                c1=c1,
                c2=c2)
  
  #Check selection constraints
  if (any(A%*%y.a>b)) {
    warning("A.1%*%y.a>b.1, conditioning is not fulfilled")
  }
  
  pvals <- rep(1,p)
  #for two sided p-values
  #pv <- 2*pmin(cdf, 1-cdf)
  #pvals[chosen] <- pv
  
  #for one-sided p-values
  for (i in 1:s){
    if (beta[chosen[i]] > 0) {
      pvals[chosen[i]] <- 1-cdf[i]
    } else {
      pvals[chosen[i]] <- cdf[i]
    }
  }

  #Give warnings if we have pvals outside of [0,1]
  out_of_range_pvals <- pvals[pvals < 0 | pvals > 1]
  if (length(out_of_range_pvals) > 0) {
    warning("Found invalid pvals:", toString(out_of_range_pvals))
    
  }
  #clip all values to [0,1], as it may happen due to numerical instabilities that some p-values fall slightly outside this range
  pvals <- pmin(pmax(pvals, 0), 1)
  
  return(list(pvals=pvals, beta_comb = beta_comb, vlo = vlo, vup = vup))
}

beta.split <- function(x, y, split, beta, sigma_squ){
  
  #Split the data
  n <- length(y)
  p <- length(beta)
  n.a <- length(split)
  n.b <- n-n.a
  x.b <- matrix(x[-split, ], nrow = n.b)
  y.b <- y[-split, ]

  chosen <-  which(abs(beta) > 0) # selected variables
  s <- length(chosen)

  if (s == 0)
    stop("0 variables were chosen by the Lasso")
  
  #extract active variables from split B
  x.Mb <- x.b[, chosen]#(n.b x s)
  
  #compute the moore penrose inverse of active variables in split B
  x.Mb.i <- ginv(x.Mb)
  
  #compute the beta_split
  beta_split<-x.Mb.i %*% y.b
  
  var_vector <- sigma_squ*diag(solve(t(x.Mb)%*%x.Mb))
  pvals <- rep(1,p)
  for (i in 1:s){
    cdf <- pnorm(beta_split[i],mean = 0, sd = sqrt(var_vector[i]))
    #for two sided tests
    pv <- 2*min(cdf, 1-cdf)
    pvals[chosen[i]] <- pv
    #for one sided tests
    # if (beta[chosen[i]] > 0) {
    #   pvals[chosen[i]] <- 1-cdf
    # } else {
    #   pvals[chosen[i]] <- cdf
    # }
  }
  return(list(pvals_split = pvals, beta_split = beta_split))
}



beta.posi <- function(x, y, split, beta, lambda,
                      sigma_squ){
  #Split the data
  n <- length(y)
  p <- length(beta)
  n.a <- length(split)
  x.a <- x[split, ]
  y.a <- y[split, ]

  #Sigma gets chosen in accordance with the distribution of y_A in Lemma 2.4.1
  Sigma <- diag(n.a)*sigma_squ
  
  chosen <-  which(abs(beta) > 0) # selected variables
  s <- length(chosen)
  b.signs <- sign(beta[chosen])
  
  if (s == 0)
    stop("0 variables were chosen by the Lasso")
  
  #extract active variables from both splits
  x.Ma <- x.a[, chosen]#(n.a x s)
  
  #extract inactive variables from both splits
  x_Ma <- x.a[, -chosen]#(n.a x (p-s))
  
  #compute the moore penrose inverse of active variables in both splits
  x.Ma.i <- ginv(x.Ma)#(s x na)
  x.Ma.ti <- ginv(t(x.Ma))
  
  #compute the projection matrix onto active columns
  p.Ma <- x.Ma %*% x.Ma.i#(n.a x n.a)
  
  beta_posi<-x.Ma.i %*% y.a
  
  #Get active affine constraints on split A
  C <- solve(crossprod(x.Ma,x.Ma))#(X_Ma^T*X_Ma)^-1
  A.1 <- -diag(x = b.signs, nrow = s) %*% x.Ma.i # (s x n.a)
  b.1 <- -lambda*diag(x = b.signs, nrow = s) %*% C %*% b.signs#Here Christoph differentiates for intercept
  
  A <- A.1
  b <- b.1
  
  #compute truncation limits
  vup <- rep(0,s)
  vlo <- rep(0,s)
  for (i in 1:s){
    eta <- x.Ma.i[i,]
    c <- (Sigma %*% eta) / as.numeric((t(eta) %*% Sigma %*% eta))
    z <- (diag(n.a) - c %*% t(eta)) %*% y.a
    den <- A%*%c
    resid <- b-A %*% z
    #We do not consider the set V^0(z) as defined in Lee et al. p.10
    ind.vup <- (den > 0)
    ind.vlo <- (den < 0)
    if (any(ind.vup)){
      vup[i] <- min(resid[ind.vup]/den[ind.vup])
    }else {
      vup[i] <- Inf
    }
    
    if (any(ind.vlo)){
      vlo[i] <- max(resid[ind.vlo]/den[ind.vlo])
    }else {
      vlo[i] <- -Inf
    }
  }
  
  tau.M <- sigma_squ
  tau.2 <- diag(tau.M*solve(t(x.Ma)%*%x.Ma))
  
  #Null hypothesis
  theta.2 <- rep(0,s)
  
  #Check selection constraints
  if (any(A%*%y.a>b)) {
    warning("A.1%*%y.a>b.1, conditioning is not fulfilled")
  }
  
  repl_for_1 = 1 - (1e-30)
  
  pvals <- rep(1,p)
  for (i in 1:s){
    cdf<-ptruncnorm(q=beta_posi[i], a=vlo[i], b=vup[i], mean=theta.2[i], sd=sqrt(tau.2[i]))

    #pnormTrunc can not handle probabilites that are very close to 1
    cdf<-ifelse(is.nan(cdf), repl_for_1, cdf)
    #for two sided p-values
    #pv <- 2*min(cdf, 1-cdf)
    #pvals[chosen[i]] <- pv
    #for one-sided p-values
    if (beta[chosen[i]] > 0) {
      pvals[chosen[i]] <- 1-cdf
    } else {
      pvals[chosen[i]] <- cdf
    }
  }
  
  return(list(pvals=pvals, beta_posi=beta_posi,vlo = vlo, vup = vup))
}