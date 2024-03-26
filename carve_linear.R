#splits the data, performs selection on one split, calculates p-values of carving estimator as in Drydale's paper

#set.seed(42)
#Set the seed to have replicabiltiy while debugging

carve.linear <- function(x, y, fraction = 0.9, args.model.selector = list(intercept = FALSE),
                         sigma=sigma){
  
  #Normalize x and y before starting:
  y<-y-mean(y)

  for (j in dim(x)[2]){
    xjbar<-mean(x[,j])
    #Calculate the variance with 1/n in the denominator as per Bühlmann's HDS lecture:
    sigma_j<-sum((x[,j]-xjbar)^2)/length(x[,j])
    for (i in dim(x)[1]){
      x[i,j]<-(x[i,j]-xjbar)/sigma_j
    }
  }
  
  #1-split-select from Christoph's carve.Lasso:
  #Splits data and selects active variables while checking constraints.
  split.select.list <- split.select(x,y,fraction = fraction)
  beta <- split.select.list$beta
  #beta <- beta[-1] #exclude intercept for now
  lambda <- split.select.list$lambda
  split <- split.select.list$split
  n <- length(y)
  p <- length(beta)
  n.a <- length(split)
  n.b <- n-n.a
  x.a <- x[split, ]
  y.a <- y[split, ]
  x.b <- x[-split, ]
  y.b <- y[-split, ]
  sigma <- sigma
  #Sigma gets chosen in accordance with the distribution of y_A in Lemma 3.2
  Sigma <- diag(n.a)*sigma
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
  
  #Check for well-definedness of moore penrose inverses, and hence also of beta_carve_D, checking if the rank of t(x.M)%*%x.M
  #is such that it allows inversion
  x.Ma.cross <- crossprod(x.Ma,x.Ma)#t(x.Ma)%*%x.Ma, s x s
  x.Mb.cross <- crossprod(x.Mb,x.Mb)#s x s
  x.Ma.cross.rank <- rankMatrix(x.Ma.cross)
  x.Mb.cross.rank <- rankMatrix(x.Mb.cross)
  if ((x.Ma.cross.rank < s) && (x.Mb.cross.rank < s)){
    cat("the dimension of t(x.Ma)%*%x.Ma is ", dim(x.Ma.cross)[1], "x",dim(x.Ma.cross)[2], "but the rank of it is ", x.Ma.cross.rank, "\n")
    cat("the dimension of t(x.Mb)%*%x.Mb is ", dim(x.Mb.cross)[1], "x",dim(x.Mb.cross)[2], "but the rank of it is ", x.Mb.cross.rank, "\n")
    stop("Both t(x.Ma)%*%X.Ma and t(x.Mb)%*%X.Mb are not invertible")
  }
  if (x.Ma.cross.rank < s){
    cat("the dimension of t(x.Ma)%*%x.Ma is ", dim(x.Ma.cross)[1], "x",dim(x.Ma.cross)[2], "but the rank of it is ", x.Ma.cross.rank, "\n")
    stop("t(x.Ma)%*%X.Ma) is not invertible")
  }
  if (x.Mb.cross.rank < s){
    cat("the dimension of t(x.Mb)%*%x.Mb is ", dim(x.Mb.cross)[1], "x",dim(x.Mb.cross)[2], "but the rank of it is ", x.Mb.cross.rank, "\n")
    stop("t(x.Mb)%*%X.Mb) is not invertible")
  }
  
  
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
  
  A <- rbind(A.0, A.1)# (2p-s) x n.a
  b <- rbind(b.0,b.1)
  

  
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
    #We do not consider the set V^0(z) as defined in Lee p.10, because
    #Drysdale does not do so either
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
      vlo[i] <- -Inf
    }
    norm_consts[i] <- v.i.norm
  }
    
  #Scale back, this is what Drysdale does, not sure if necessary
  eta_var <- sigma * (norm_consts^2)
  vlo <- vlo * norm_consts
  vup <- vup * norm_consts
  # Turn all signs positive, as Drysdale does
  mask = (b.signs == -1)
  vlo[mask] <- -vlo[mask]
  vup[mask] <- -vup[mask]
  #See comments in the RMD file for tau.1 and tau.2
  #tau.1 <- sigma
  #tau.2 <- eta_var
  

  #REMARK: Drysdale sets theta1 = theta2 = beta_null for the sntn dist, where beta_null is the assumed beta under the null, 
  #so in our case an all zeros vector of dimension beta_carve_D, for reference: see parameters of run_inference in _lasso.py
  theta.1 <- rep(0,s)
  theta.2 <- theta.1
  
  #y~N(x beta^0, tau^2 I_n)
  tau.M=sigma

  #Defined beta^M=0 for testing the null hypothesis
  beta.M=rep(0,s)  


  pvals<-1-SNTN_CDF(z=beta_carve_D,
           mu1=theta.1,
           tau1=diag(tau.M*solve(t(x.Mb)%*%x.Mb)),
           mu2=theta.2,
           tau2=diag(tau.M*solve(t(x.Ma)%*%x.Ma)),
           a=vlo,
           b=vup,
           c1=c1,
           c2=c2)
  
  
  return(list(pvals = pvals,split = split, beta = beta, lambda = lambda))
}



#----------------------For manually computing these values: -----------------------------------------
#In the Drysdale paper, sigma_1 gets defined via the entry jj for all j. Therefore we take the diagonal
#of the matrix and work with that
#Since it's a diagonal matrix, we can just apply the inverse at the end, which is
#computationally easier
# sigma.1=(tau.M/n^2)*((n.b^2*diag((t(x.Mb)%*%x.Mb))^-1)+(n.a^2*diag((t(x.Ma)%*%x.Ma))^-1))
# sigma.2=tau.M*diag((t(x.Ma)%*%x.Ma))^-1
# w=(vlo-beta.M)/(sqrt(tau.M)*diag((t(x.Ma)%*%x.Ma))^(-1/2))
# delta=(vup-beta.M)/(sqrt(tau.M)*diag((t(x.Ma)%*%x.Ma))^(-1/2))
#I'm assuming that we need the jj-th entry again here everywhere. Otherwise division 
# wouldn't make sense and we'd need to build some sort of framework for allowing the sqrt
# rho=sqrt(n*n.a)*(diag(t(x.Ma)%*%x.Ma))^(-1/2)/
#   sqrt(n.b^2*diag((t(x.Mb)%*%x.Mb))^-1 + n.a^2*diag((t(x.Ma)%*%x.Ma))^-1)
