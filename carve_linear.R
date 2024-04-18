
#' splits the data, calculates p-values of carving estimator as in Drydale's paper
#'
#' @param x (matrix dim: n x p) the full design matrix
#' @param y (vector dim: n) the full response vector
#' @param split (vector dim: n_a) indices of observations used for selection
#' @param beta (vector dim: p) full beta obtained from lasso selection
#' @param lambda (numeric) lambda from lasso selection
#' @param sigma (numeric) true variance of data generating process y~N(0,sigma^2*I_n)
#' @param normalize_truncation_limits (bool) if true, truncation limits are normalized
#'
#' @return list containing p-values of Drysdales carving estimator and all of the inputs used for SNTN_distribution, as well
#' as the norms of the directions of interest during the computation of the truncation limits

carve.linear <- function(x, y, split, beta, lambda,
                         sigma=sigma, normalize_truncation_limits = FALSE){
  

  #Split the data
  n <- length(y)
  p <- length(beta)
  n.a <- length(split)
  n.b <- n-n.a
  x.a <- x[split, ]
  y.a <- y[split, ]
  x.b <- x[-split, ]
  y.b <- y[-split, ]
  sigma_squ <- sigma
  
  #Sigma gets chosen in accordance with the distribution of y_A in Lemma 3.2
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
  beta_carve_D <- ((n.a/n)*x.Ma.i %*% y.a + (n.b/n)*x.Mb.i %*% y.b)
  
  #Get inactive affine constraints on split A, as this is the group where we are doing POSI(Lee et al. page 8)
  A.0.up <- (t(x_Ma) %*% (diag(n.a) - p.Ma))#(p-s) x n.a
  A.0 <- 1/lambda*rbind(A.0.up, -A.0.up) #2*(p-s) x n.a
  b.0.temp <-t(x_Ma) %*% x.Ma.ti %*% b.signs
  b.0.up <- rep(1,p-s) - b.0.temp
  b.0.lo <- rep(1,p-s) + b.0.temp
  b.0 <- rbind(b.0.up, b.0.lo)#2*(p-s) x n.a
  
  #Get active affine constraints on split A
  C <- solve(crossprod(x.Ma,x.Ma))#(X_Ma^T*X_Ma)^-1
  A.1 <- -diag(x = b.signs, nrow = s) %*% x.Ma.i # (s x n.a)
  b.1 <- -lambda*diag(x = b.signs, nrow = s) %*% C %*% b.signs#Here Christoph differentiates for intercept
  
  #A <- rbind(A.0, A.1)# (2p-s) x n.a
  #b <- rbind(b.0,b.1)
  A <- A.1
  b <- b.1
  
  #Following a mix of Lee (https://github.com/selective-inference/R-software/blob/master/selectiveInference/R/funs.fixed.R) from line 231
  #and Drysdale (https://github.com/ErikinBC/sntn/blob/main/sntn/_lasso.py) from line 195
  vup <- rep(0,s)
  vlo <- rep(0,s)
  norm_consts <- rep(0,s)
  for (i in 1:s){
    v.i <- x.Ma.i[i,]
    v.i.norm <- sqrt(sum(v.i^2))
    if(normalize_truncation_limits){
      eta <- b.signs[i]*v.i/v.i.norm
    }
    else {
      eta <- x.Ma.i[i,]
    }
    c <- (Sigma %*% eta) / as.numeric((t(eta) %*% Sigma %*% eta))
    z <- (diag(n.a) - c %*% t(eta)) %*% y.a
    den <- A%*%c
    resid <- b-A %*% z
    #We do not consider the set V^0(z) as defined in Lee p.10, because
    #Drysdale does not do so either
    ind.vup <- (den > 0)
    #The next line ensures that when calculating vup we only consider the cases where resid is positive
    #ind.vup <- which(ind.vup == TRUE)[which(resid[ind.vup]>0)]
    ind.vlo <- (den < 0)
    #Same here for vlo, this is important for vlo to be smaller than vup
    #ind.vlo <- which(ind.vlo == TRUE)[which(resid[ind.vlo]>0)]
    ind.v0 <- (den == 0)
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
    norm_consts[i] <- v.i.norm
  }
  
  tau.M <- sigma_squ
  neg_mask = (b.signs == -1)
  
  if(normalize_truncation_limits){
    #Scale back, this is what Drysdale does, not sure if necessary
    eta_var <- sigma_squ * (norm_consts^2)
    vlo <- vlo * norm_consts
    vup <- vup * norm_consts
    
    # Turn all signs positive, as Drysdale does
    #NEW CHANGES: Drysdales command V[mask] = -V[mask][:,[1,0]] also swaps vlo and vup at the positions of mask, not only changing their signs
    vlo_temp <- vlo
    vlo[neg_mask] <- -vup[neg_mask]
    vup[neg_mask] <- -vlo_temp[neg_mask]
  }
  else{
    eta_var <- tau.M
  }
  #See comments in the RMD file for tau.1 and tau.2
  tau.1 <- diag(tau.M*solve(t(x.Mb)%*%x.Mb))
  tau.2 <- diag(eta_var*solve(t(x.Ma)%*%x.Ma))
  
  
  #REMARK: Drysdale sets theta1 = theta2 = beta_null for the sntn dist, where beta_null is the assumed beta under the null, 
  #so in our case an all zeros vector of dimension beta_carve_D, for reference: see parameters of run_inference in _lasso.py
  theta.1 <- rep(0,s)
  theta.2 <- rep(0,s)
  
  cdf<-SNTN_CDF(z=beta_carve_D,
                mu1=theta.1,
                tau1=tau.1,
                mu2=theta.2,
                tau2=tau.2,
                a=vlo,
                b=vup,
                c1=c1,
                c2=c2)
  
  #Inserted this to check for direction of pvals, if the sign of beta_select is negative, we take cdf
  #to be the pv, else we take 1-cdf(line 319 in _lasso.py from Drysdales code)
  #pv <- ifelse(neg_mask, cdf, 1-cdf)
  #For a two sided test we could also take this, as i dont fully understand the above
  pv <- 2*pmin(cdf, 1-cdf)
  pvals <- rep(1,p)
  pvals[chosen] <- pv
  
  #Give warnings if we have pvals outside of [0,1]
  out_of_range_pvals <- pvals[pvals < 0 | pvals > 1]
  if (length(out_of_range_pvals) > 0) {
    warning("Found invalid pvals:", toString(out_of_range_pvals))
    
  }
  #clip all values to [0,1]
  pvals <- pmin(pmax(pvals, 0), 1)
  
  return(list(pvals=pvals, norm_consts=norm_consts, beta_carve_D = beta_carve_D, tau.1 = tau.1,
              tau.2 = tau.2, vlo = vlo, vup = vup, c1=c1, c2=c2, x.Ma.i=x.Ma.i, x.Mb.i=x.Mb.i,
              y.a=y.a, y.b=y.b))
}