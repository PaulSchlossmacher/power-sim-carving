
SNTN_CDF <- function(z,mu1, tau1, mu2, tau2, a, b, c1, c2) {
  #Following Lemma 3.1 from Drysdale paper
  theta1 <- c1 * mu1 + c2 * mu2
  sigma1 <- (c1)^2 * (tau1)^2 + (c2)^2 * (tau2)^2
  r.sigma1 <- sqrt(sigma1)
  theta2 <- mu2
  sigma2 <- tau2
  r.sigma2 <- sqrt(sigma2)
  rho <- c2 * r.sigma2 / r.sigma1
  if (rho <= -1 | rho >= 1) {
    stop("Invalid correlation coefficient rho. It should be between -1 and 1.")
  }
  lambda <- rho / sqrt(1 - rho^2)
  gamma <- lambda/rho
  w <- (a-theta2)/r.sigma2
  delta <- (b-theta2)/r.sigma2
  m1_z <- (z-theta1)/r.sigma1
  
  #This is the covariance matrix for the bivariate normal dist used by Drysdale
  cov_matrix <- matrix(c(sigma1, rho * r.sigma1*r.sigma2,
                    rho * r.sigma1*r.sigma2, sigma2), nrow = 2)
  mean_delta <- c(m1_z,delta)
  mean_w <- c(m1_z,w)
  #lower and upper are integration bounds when calculating bivariate CDF
  bn_cdf_delta <- pmvnorm(lower = c(-Inf, -Inf), upper = mean_delta, mean = c(0,0), sigma = cov_matrix,keepAttr = FALSE)
  bn_cdf_w <- pmvnorm(lower = c(-Inf, -Inf), upper = mean_w, mean = c(0,0), sigma = cov_matrix,keepAttr = FALSE)
  n_cdf_delta <- pnorm(delta)
  n_cdf_w <- pnorm(w)
  sntn_cdf <- (bn_cdf_delta - bn_cdf_w)/(n_cdf_delta - n_cdf_w)
  return(sntn_cdf)
}


