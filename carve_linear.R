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
  Sigma <- diag(n.a)*sigma
  
  chosen <-  which(abs(beta) > 0) # selected variables
  s <- length(chosen)
  b.signs <- sign(beta[chosen])
  if (s == 0) {
    return(NULL)
  }
  
  
  #extract active variables from both splits
  x.Ma <- x.a[, chosen]
  x.Mb <- x.b[, chosen]
  
  #extract inactive variables from both splits
  x_Ma <- x.a[, -chosen]#(n.a x (p-s))
  x_Mb <- x.b[, -chosen]
  
  
  #compute the moore penrose inverse of active variables in both splits
  x.Ma.i <- ginv(x.Ma)
  x.Ma.ti <- ginv(t(x.Ma))
  x.Mb.i <- ginv(x.Mb)
  
  #compute the projection matrix onto active columns
  p.Ma <- x.Ma %*% x.Ma.i#(n.a x n.a)
  
  #compute the beta_carve from drysdales paper
  beta_carve_D <- n*(fraq*x_Ma.i %*% y_a + (1-fraq)*x_Mb.i %*% y_b)
  
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
  c <- Sigma %*% t(x.Ma.i) %*% solve(x.Ma.i %*% Sigma %*% t(x.Ma.i))
  
  z <- (diag(n.a) - c %*% x.Ma.i) %*% y.a
  
  ind.vup <- (A %*% c > 0)
  ind.vlo <- (A %*% c < 0)
  
  
  vup <- max((b-A %*% z)/(A %*% c)[ind.vup])
  vlo <- max((b-A %*% z)/(A %*% c)[ind.vlo])
  
  
}
