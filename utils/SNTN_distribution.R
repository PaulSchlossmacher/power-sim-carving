
#' SNTN cumulative distribution function, following Lemma 2.4.1
#'
#' @param z (vector): values at which the sntn_cdf is evaluated
#' @param mu1 (vector): means of the normal cdf for each entry of z
#' @param tau1 (vector): variances of the normal cdf for each entry of z
#' @param mu2 (vector): means of the truncated normal cdf for each entry of z
#' @param tau2 (vector): variances of the truncated normal cdf for each entry of z
#' @param a (vector): lower truncation limits for the truncated normal cdf 
#' @param b (vector): upper truncation limits for the truncated normal cdf
#' @param c1 (numeric): coefficient indicating the influence of the normal cdf in the sntn cdf
#' @param c2 (numeric): coefficient indicating the influence of the truncated normal cdf in the sntn cdf
#'
#' @return P(x<z) under the sntn cdf P

SNTN_CDF <- function(z,mu1, tau1, mu2, tau2, a, b, c1, c2) {
  
  if (!all(a<b)){
    warning("We encountered ", toString(sum(!(a<b))), "/", toString(length(z)), " entries of vlo that are larger than vup")
  }
  s <- length(z)
  theta1 <- c1 * mu1 + c2 * mu2
  sigma1 <- (c1)^2 * tau1 + (c2)^2 * tau2
  r.sigma1 <- sqrt(sigma1)
  theta2 <- mu2
  sigma2 <- tau2
  r.sigma2 <- sqrt(sigma2)
  rho <- c2 * r.sigma2 / r.sigma1
  if (min(rho) < -1 | max(rho) > 1) {
    stop("Invalid correlation coefficient rho. It should be between -1 and 1.")
  }
  lambda <- rho / sqrt(1 - rho^2)
  gamma <- lambda/rho
  w <- (a-theta2)/r.sigma2
  delta <- (b-theta2)/r.sigma2
  m1_z <- (z-theta1)/r.sigma1
  
  sntn_cdf_arr <- rep(0,s)
  for (i in 1:s){
    #Standard bivariate normal with correlation rho:
    cov_matrix <- matrix(c(1, rho[i],
                           rho[i], 1), nrow = 2)
    mean_delta <- c(m1_z[i],delta[i])
    mean_w <- c(m1_z[i],w[i])
    #lower and upper are integration bounds when calculating bivariate CDF
    bn_cdf_delta <- pmvnorm(lower = c(-Inf, -Inf), upper = mean_delta, mean = c(0,0), sigma = cov_matrix,keepAttr = FALSE)
    bn_cdf_w <- pmvnorm(lower = c(-Inf, -Inf), upper = mean_w, mean = c(0,0), sigma = cov_matrix,keepAttr = FALSE)
    n_cdf_delta <- pnorm(delta[i])
    n_cdf_w <- pnorm(w[i])

    #If clauses to handle division by 0: This is slightly different to what Drysdale
    #does, but the most conservative approach (in the sense of avoiding Type I errors)
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

#This is from Lemma 3.1 in Drysdales paper
SNTN_PDF <- function(z,mu1, tau1, mu2, tau2, a, b, c1, c2) {
  s <- length(z)
  theta1 <- c1 * mu1 + c2 * mu2
  sigma1 <- (c1)^2 * tau1 + (c2)^2 * tau2
  r.sigma1 <- sqrt(sigma1)
  theta2 <- mu2
  sigma2 <- tau2
  r.sigma2 <- sqrt(sigma2)
  rho <- c2 * r.sigma2 / r.sigma1
  if (min(rho) < -1 | max(rho) > 1) {
    stop("Invalid correlation coefficient rho. It should be between -1 and 1.")
  }
  lambda <- rho / sqrt(1 - rho^2)
  gamma <- lambda/rho
  w <- (a-theta2)/r.sigma2
  delta <- (b-theta2)/r.sigma2
  m1_z <- (z-theta1)/r.sigma1
  
  sntn_pdf_arr <- rep(0,s)
  for (i in 1:s){
    num <- dnorm(m1_z[i])*(pnorm((gamma*delta - lambda*m1_z)[i]) - pnorm((gamma*w - lambda*m1_z)[i]))
    den <- r.sigma1[i]*(pnorm(delta[i])-pnorm(w[i]))
    sntn_pdf <- num/den
    sntn_pdf_arr[i] <- sntn_pdf
    
  }
  return(sntn_pdf_arr)
}
