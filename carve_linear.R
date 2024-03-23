#splits the data, performs selection on one split, calculates p-values of carving estimator as in Drydale's paper




carve.linear <- function(x, y, fraction = 0.9, FWER = TRUE, family = "gaussian", model.selector = lasso.cvcoef,
                          args.model.selector = list(intercept = TRUE, standardize = FALSE, tol.beta = 1e-5),
                          df.corr = FALSE, args.lasso.inference = list(sigma = sigma), verbose = FALSE){
  #binomial not implemented yet
  if (!(family %in% c("gaussian", "binomial")))
    stop ("Invalid family provided, can only deal with gaussian and binomial")
  
  args.model.selector$family <- family
  args.lasso.inference$family <- family
  
  split.select.list <- split.select(x,y,fraction = fraq)
  beta <- split.select.list$beta
  beta <- beta[-1]#exclude intercept for now
  lambda <- split.select.list$lambda
  split <- split.select.list$split
  n <- length(y)
  p <- length(beta)
  n.a <- length(split)
  n.b <- n-n.a
  x.a <- x[split, ]
  y.a <- y[split,]
  x.b <- x[-split, ]
  y.b <- y[-split]
  sigma <- args.lasso.inference$sigma
  Sigma <- diag(n.a)*sigma#cov of y
  c1 <- 1-fraq
  c2 <- fraq
  
  chosen <-  which(abs(beta) > 0) # selected variables
  s <- length(chosen)
  b.signs <- sign(beta[chosen])
  if (s == 0) {
    return(NULL)
  }
  
  #extract active variables from both splits
  x.Ma <- x.a[, chosen]#(n.a x s)
  x.Mb <- x.b[, chosen]
  
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
  beta_carve_D <- n*((n.a/n)*x.Ma.i %*% y.a + (n.b/n)*x.Mb.i %*% y.b)
  
  #Get inactive affine constraints on split A, as this is the group where we are doing POSI(Lee et al. page 8)
  A.0.up <- (t(x_Ma) %*% (diag(n.a) - p.Ma))#(p-s) x n.a
  A.0 <- 1/lambda*rbind(A.0.up, -A.0.up) #2*(p-s) x n.a
  b.0.temp <-t(x_Ma) %*% x.Ma.ti %*% b.signs
  b.0.up <- rep(1,p-s) - b.0.temp
  b.0.lo <- rep(1,p-s) + b.0.temp
  b.0 <- rbind(b.0.up, b.0.lo)#2*(p-s) x n.a
  
  #Get active affine costraints on split A
  C <- solve(crossprod(x.Ma,x.Ma))#(X_Ma^T*X_Ma)^-1
  A.1 <- -diag(x = b.signs, nrow = s) %*% x.Ma.i # (s x n.a)
  b.1 <- -diag(x = b.signs, nrow = s) %*% C %*% (lambda * b.signs)#Here Christoph differentiates for intercept
  
  A <- rbind(A.0, A.1)# (2p-s) x n.a
  b <- rbind(b.0,b.1)
  #c <- Sigma %*% t(x.Ma.i) %*% solve(x.Ma.i %*% Sigma %*% t(x.Ma.i))
  #z <- (diag(n.a) - c %*% x.Ma.i) %*% y.a
  
  #TODO: Need to handle these if empty, set to inf or -inf respectively
  # ind.vup <- (A %*% c > 0)
  # ind.vlo <- (A %*% c < 0)
  
  #TODO:There's an error here, as the dimensions do not match. I think we have overseen that A%*%c is a matrix, hence division is ambiguous.
  #Drysdale solves this issue over a for loop in inference_on_screened inside _lasso.py. There he does it row by row of eta.T.
  #If i am not mistaken, then v_i is a column of c (equation 5.3 in Lee et al.) i am just not sure why he uses the signs of beta_hat in this equation, 
  #while resid_i is z.
  # vup <- min((b-A %*% z)/(A %*% c)[ind.vup])
  # vlo <- max((b-A %*% z)/(A %*% c)[ind.vlo])
  
  
  #Following a mix of Lee (https://github.com/selective-inference/R-software/blob/master/selectiveInference/R/funs.fixed.R) from line 231
  #and Drysdale (https://github.com/ErikinBC/sntn/blob/main/sntn/_lasso.py) from line 195
  vup <- rep(0,s)
  vlo <- rep(0,s)
  norm_consts <- rep(0,s)
  for (i in 1:s){
    v.i <- x.Ma.i[i,]
    v.i.norm <- sqrt(sum(v.i^2))
    eta <- b.signs[i]*v.i/v.i.norm
    c <- (Sigma %*% eta) / as.numeric((t(eta) %*% Sigma %*% eta))
    z <- (diag(n.a) - c %*% t(eta)) %*% y.a
    den <- A%*%c
    resid <- b-A %*% z
    ind.vup <- (A %*% c > 0)
    ind.vlo <- (A %*% c < 0)
    if (any(ind.vup)){
      vup[i] <- min(resid[ind.vup]/den[ind.vup])
    }else {
      vup[i] <- Inf
    }
    
    if (any(ind.vlo)){
      vlo[i] <- max(resid[ind.vlo]/den[ind.vlo])
    }else {
      vlo[i] <- Inf
    }
    norm_consts[i] <- v.i.norm
  }
    
  #Scale back, this is what Drysdale does, not sure if necessary
  eta_var <- sigma * (norm_consts^2)
  vlo <- vlo * norm_consts
  vup <- vup * norm_consts
  # Return to original signs
  mask = (b.signs == -1)
  vlo[mask] <- -vlo[mask]
  vup[mask] <- -vup[mask]
  tau.1 <- sigma
  tau.2 <- eta_var
  
  
  
  #I'm not sure about this yet, but at least in the analogy to OLS I think it makes sense:
  #There we also build confidence intervals etc. based on beta^hat.
  beta.M=beta_carve_D
  #REMARK: Drysdale sets theta1 = theta2 = beta_null for the sntn dist, where beta_null is the assumed beta under the null, 
  #so in our case an all zeros vector of dimension beta_carve_D, for reference: see parameters of run_inference in _lasso.py
  theta.1 <- rep(0,length(beta_carve_D))
  theta.2 <- theta.1
  ################## TEST SNTN_distribution################################
  mu1 <- rep(0, length(beta_carve_D))
  mu2 <-  mu1
  a <- vlo
  b <- vup
  sntn_cdf_arr <- SNTN_CDF(beta_carve_D,mu1, tau.1, mu2, tau.2, a, b, c1, c2)
  #Liefert alles 1, da mean_delta und mean_w riesig sind, vielleicht müssten wir irgendwo noch unsere daten normalisieren, 
  #Weiss nicht ob das so sein soll, sonst hätten wir alle pvalues gleich 0
  #############################################################################
  
  #y~N(x beta^0, tau^2 I_n)
  tau.M=sigma
  #We assume - to start off with - that tau^2=tau_M^2. Drysdale does mention on p. 3 that 
  #it could happen that sigma_M^2>sigma^2 due to some covariates being left out
  #REMARK: Drysdale chooses our tau.M for tau.1(for the distribution of beta_split), but uses some scaled version for 
  #the tau.2 (for the truncated distribution of beta_posi). The choice of this scaling is not yet clear to me
  #For refernce: see lasso.run_inference in run_carve, in the if clause he defines tau22 = eta2var, which is a scaled version of 
  #sigma2
  

  #In the Drysdale paper, sigma_1 gets defined via the entry jj for all j. Therefore we take the diagonal
  #of the matrix and work with that
  #Since it's a diagonal matrix, we can just apply the inverse at the end, which is
  #computationally easier
  sigma.1=(tau.M/n^2)*((n.b^2*diag((t(x.Mb)%*%x.Mb))^-1)+(n.a^2*diag((t(x.Ma)%*%x.Ma))^-1))
  sigma.2=tau.M*diag((t(x.Ma)%*%x.Ma))^-1
  w=(vlo-beta.M)/(sqrt(tau.M)*diag((t(x.Ma)%*%x.Ma))^(-1/2))
  delta=(vup-beta.M)/(sqrt(tau.M)*diag((t(x.Ma)%*%x.Ma))^(-1/2))
  
  #I'm assuming that we need the jj-th entry again here everywhere. Otherwise division 
  #wouldn't make sense and we'd need to build some sort of framework for allowing the sqrt
  rho=sqrt(n*n.a)*(diag(t(x.Ma)%*%x.Ma))^(-1/2)/
    sqrt(n.b^2*diag((t(x.Mb)%*%x.Mb))^-1 + n.a^2*diag((t(x.Ma)%*%x.Ma))^-1)
  #This return is used for debugging
  return(list(sigma.1 = sigma.1, sigma.2 = sigma.2, w = w, delta = delta, w = w))
  #REMARK: For my definition of sntn_cdf we dont need the explicit sigma.1,sigma.2, w, delta, rho, as they get calculated above.
  #It seems to me, that my sntn_cdf would deliver different results for these quantities, see e.g. sigma2 <- tau2 in sntn_cdf, 
  #wheras sigma.2 in the lines above takes into account the whole variance of beta_posi. A general seperate sntn_cdf function needs a 
  #bit more refinement to perform as desired. So we can try to do it inside of here again:

  # pvals <- 1 - SNTN_CDF(beta_carve_D, ...)
  
}



### Theoretical notes:

#I asked Filip on Friday how you actually compute things when you only want to condition on one sign pattern.
#Lee makes this clear on p. 15:
#"Conditioning on the signs means that we only have to compute the interval [V−s (z), V+s (z)]
#for the sign pattern s that was actually observed."

#We see right under Theorem 5.3 in Lee, that V-s(z) and V+s(z) are defined through
# A=As and b=bs
#And s influences the definitions of A1(M,s) and b1(M,s) respectively.

#Since s is in {-1,1}^|M|, it's only defined for variables that are actually selected, so
#the computation of the signs is straightforward (I mention this, because we had some confusion
#with a similar thing in another paper where we had s in [-1,1]^|M| or sth like this )


#Question: Which beta^hat are we actually using though to get the signs?
#A priori all of beta^carve, beta^POSI and beta^SPLIT seem at least viable

#Thinking about it, I guess that since we are talking about M (i.e. M_A) all the time, it is probably 
#beta^Split, which is also the beta we are working with in the code above. In fact, Filip already 
#implemented it exactly like that above.


#Question: If we only have eta in R^nx1 for a single polyhedron and eta_M in R^nx|M| for the union of polyhedra:
#What eta_M do we actually use now when we additionally condition on the signs, to only have one polyhedron?

3+5