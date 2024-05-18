rm(list = ls())


#Local, user specific path, that should work for both of us:
Local_path<-getwd()
hdi_adjustments_path<-paste(Local_path, "/Multicarving-Christoph/inference/hdi_adjustments.R", sep="")
carving_path<-paste(Local_path, "/Multicarving-Christoph/inference/carving.R", sep="")
sample_from_truncated_path<-paste(Local_path, "/Multicarving-Christoph/inference/sample_from_truncated.R", sep="")
tryCatchWE_path<-paste(Local_path, "/Multicarving-Christoph/inference/tryCatch-W-E.R", sep="")

#Different paths here, because they're "our own" functions
SNTN_distribution_path<-paste(Local_path, "/SNTN_distribution.R", sep="")
split_select_function_path<-paste(Local_path, "/split_select.R", sep="")
carve_linear_path<-paste(Local_path, "/carve_linear.R", sep="")


library(MASS)
library(mvtnorm)
library(glmnet)
library(Matrix)
library(tictoc)
library(hdi)
library(selectiveInference)
library(doSNOW)
library(parallel)
library(doRNG)
library(truncnorm)
library(git2r)
library(ggplot2)


source(hdi_adjustments_path)
source(carving_path)
source(sample_from_truncated_path)
source(tryCatchWE_path)
source(SNTN_distribution_path)
source(split_select_function_path)
source(carve_linear_path)


#Short general test,taking N/TN parameters to be same with large truncation limits gives a pretty normal cdf as expected
#To show that it can be much different, I increased the variance of the TN and gave it small truncation limits. This way we get a very wide gaussian,
#which falls fast at the truncation limits. To see the effect of the TN, i added the contribution of the normal cdf to the SNTN in red
test_points <- seq(-4,4, 0.01)
theta1 <- rep(0,length(test_points))
theta2 <- rep(0,length(test_points))
sigma1 <- rep(1,length(test_points))
sigma2 <- rep(10,length(test_points))
c1 <- 0.3
c2 <- 0.7
a <- -3
b <- 3
normal_pdf <- dnorm(test_points, mean = 0, sd = 1)
cdf<-SNTN_CDF(z=test_points,
              mu1=theta1,
              tau1=sigma1,
              mu2=theta2,
              tau2=sigma2,
              a=a,
              b=b,
              c1=c1,
              c2=c2)
pdf<-SNTN_PDF(z=test_points,
              mu1=theta1,
              tau1=sigma1,
              mu2=theta2,
              tau2=sigma2,
              a=a,
              b=b,
              c1=c1,
              c2=c2)
par(mfrow = c(1,2))
plot(test_points,pdf, main = paste0("SNTN pdf with parameters of normal distr (", theta1[1],",", sigma1[1],") and
                        of truncated normal (",theta2[1],",",sigma2[1],",",a,",",b,")" ))
lines(test_points, 0.3*normal_pdf, type = "l", col = "red")
plot(cdf, main = "SNTN cdf")




######################## Test with our relevant data #######################################
n <- 100
p <- 200
rho <- 0.6
fraq = 0.7
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
#More active variables than observations in Group B after the split:
act.index <- c(1, 5, 10, 15, 20)#active predictors
beta <- rep(0, p)
beta[act.index] <- 1
sparsity <- 5
set.seed(42) 
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta
SNR <- 1.713766 # value created for Toeplitz 0.6
sigma <- 1.2

counter <- 0
select.again <- TRUE
select.again.counter = 0
while(select.again){
  if (select.again.counter > 50){
    stop("Tried to many selection events and not one of them was conformable for beta_Drysdale")
  }
  select.again <- FALSE
  set.seed(counter)
  counter <- counter + 1
  y <- y.true + sigma * rnorm(n)
  split.select.list <- split.select(x,y,fraction = fraq)
  beta_tmp <- split.select.list$beta
  
  if(sum(beta_tmp!=0)==0){
    select.again <- TRUE
    print("0 variables where chosen by the lasso, repeating selection")
  }
  #cat("We selected ", sum(beta_tmp!=0), " predictors.\n")
  lambda <- split.select.list$lambda
  split <- split.select.list$split
  if(sum(beta_tmp!=0)>min(n*fraq, n*(1-fraq))){
    select.again <- TRUE
    select.again.counter <- select.again.counter + 1
    print("Need to split again because we selected more variables than beta_D can handle")
  }
  
}
test_points <- seq(-2,2, 0.01)
carve_D <- carve.linear(x,y,split = split, beta = beta_tmp, lambda = lambda,
                        sigma=sigma, normalize_truncation_limits = TRUE)
print(carve_D)

beta_carve_D <- carve_D$beta_carve_D
tau.1 <- carve_D$tau.1
tau.2 <- carve_D$tau.2
vlo <- carve_D$vlo
vup <- carve_D$vup
c1 <- carve_D$c1
c2 <- carve_D$c2

#Look at this:
print(vlo<vup)

