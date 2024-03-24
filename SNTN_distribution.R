
SNTN_CDF <- function(z,mu1, tau1, mu2, tau2, a, b, c1, c2) {
  #Following Lemma 3.1 from Drysdale paper
  #REMARK: Adapted to work for vector inputs, in this case it will also return a
  #vector of cdf values for each entry of z
  s <- length(z)
  theta1 <- c1 * mu1 + c2 * mu2
  #Attention here Filip wrote (tau1)^2 and (tau2)^2, but they are already passed as squares to the function
  sigma1 <- (c1)^2 * tau1 + (c2)^2 * tau2
  r.sigma1 <- sqrt(sigma1)
  theta2 <- mu2
  sigma2 <- tau2
  r.sigma2 <- sqrt(sigma2)
  rho <- c2 * r.sigma2 / r.sigma1
  if (min(rho) <= -1 | max(rho) >= 1) {
    stop("Invalid correlation coefficient rho. It should be between -1 and 1.")
  }
  lambda <- rho / sqrt(1 - rho^2)
  gamma <- lambda/rho
  w <- (a-theta2)/r.sigma2
  delta <- (b-theta2)/r.sigma2
  m1_z <- (z-theta1)/r.sigma1
  
  #This is the covariance matrix for the bivariate normal dist used by Drysdale
  #https://github.com/ErikinBC/sntn/blob/main/sntn/_cdf_bvn/_utils.py
  sntn_cdf_arr <- rep(0,s)
  for (i in 1:s){
    cov_matrix <- matrix(c(1, rho[i],
                           rho[i], 1), nrow = 2)
    mean_delta <- c(m1_z[i],delta[i])
    mean_w <- c(m1_z[i],w[i])
    #lower and upper are integration bounds when calculating bivariate CDF
    bn_cdf_delta <- pmvnorm(lower = c(-Inf, -Inf), upper = mean_delta, mean = c(0,0), sigma = cov_matrix,keepAttr = FALSE)
    bn_cdf_w <- pmvnorm(lower = c(-Inf, -Inf), upper = mean_w, mean = c(0,0), sigma = cov_matrix,keepAttr = FALSE)
    n_cdf_delta <- pnorm(delta[i])
    n_cdf_w <- pnorm(w[i])
    #TODO: Need to handle division by zero, on some runs i get all bn_cdfs and n_cdfs to be 1, which yields Nan values for sntn_cdf
    sntn_cdf <- (bn_cdf_delta - bn_cdf_w)/(n_cdf_delta - n_cdf_w)
    sntn_cdf_arr[i] <- sntn_cdf
    #TODO: Need to handle division by zero, on some runs I get all bn_cdfs and n_cdfs to be 1, which yields Nan values for sntn_cdf
    #Paul: Done, see below:
    if (bn_cdf_delta==bn_cdf_w & n_cdf_delta==n_cdf_w){
      sntn_cdf_arr[i] <- 0
    }
    else if (n_cdf_delta==n_cdf_w){
      sntn_cdf_arr[i] <- 1  
    }
    else if (bn_cdf_delta==bn_cdf_w){
      sntn_cdf_arr[i] <- 0
    }
    else{
      sntn_cdf <- (bn_cdf_delta - bn_cdf_w)/(n_cdf_delta - n_cdf_w)
      sntn_cdf_arr[i] <- sntn_cdf  
    }
    
  }
  return(sntn_cdf_arr)
}



#-------------------       Testing for the if else statements that avoid division by 0. ----------------------------------
# Res<-carve.linear(x,y,fraq, args.lasso.inference = args.lasso.inference)
# delta<-Res$delta
# w<-Res$w
# sigma1<-Res$sigma.1
# theta1<-0
# m1_z<-(Res$beta-0)/sigma1
# rho=Res$rho
# cov_matrix <- matrix(c(1, rho[1],
#                        rho[1], 1), nrow = 2)
# #Test out the denominator of F, to see whether it'S actually sensible to set the result to 0 or whether the numerator is also 0:
# mean_delta <- c(m1_z[1],delta[1])
# mean_w <- c(m1_z[1],w[1])
# #lower and upper are integration bounds when calculating bivariate CDF
# bn_cdf_delta <- pmvnorm(lower = c(-Inf, -Inf), upper = mean_delta, mean = c(0,0), sigma = cov_matrix,keepAttr = FALSE)
# bn_cdf_w <- pmvnorm(lower = c(-Inf, -Inf), upper = mean_w, mean = c(0,0), sigma = cov_matrix,keepAttr = FALSE)
# #Turns out that both summands of the numerator=1, so the numerator is also =0.
# #Therefore the more conservative decision is to set the prob to 0, so that the p-value is 1.



# #What result does the cdf give for very "average" inputs?
# #We'd hope for sth like 1/2 here.
# SNTN_CDF(z=0,
#          mu1=0,
#          tau1=1,
#          mu2=0,
#          tau2=1,
#          a=-5,
#          b=5,
#          c1=1/10,
#          c2=9/10)
# 
# #This does indeed yield 1/2, which is promising.
# 
# 
# SNTN_CDF(z=1,
#          mu1=0,
#          tau1=1,
#          mu2=0,
#          tau2=1,
#          a=-5,
#          b=5,
#          c1=1/10,
#          c2=9/10)
# #For z=1 we get 0.865, which intuitively seems quite high, but not completely unreasonable
# 
# SNTN_CDF(z=-1,
#          mu1=0,
#          tau1=1,
#          mu2=0,
#          tau2=1,
#          a=-5,
#          b=5,
#          c1=1/10,
#          c2=9/10)
# #For z=1 we get 0.134, which pairs nicely with the result for z=1, seems like symmetry around 0 is fulfilled.