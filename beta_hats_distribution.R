#Check the distribution of both beta_split and beta_posi

#Trying to test the distribution of pvalues under screening
# Clear all variables
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



#-------------------- Toeplitz Carving simulation from Christoph ----------------------
n <- 500
p <- 600
rho <- 0.6
fraq = 0.6
#fraq.vec <- c(0.7)
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
#More active variables than observations in Group B after the split:
act.index <- round(seq(from=1, to=p, length.out=70),0)#active predictors
beta <- rep(0, p)
beta[act.index] <- 1
set.seed(41)
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta
sigma_squ <- 1.2
y <- y.true + sigma_squ * rnorm(n)

#Normalize y:
y<-(y-mean(y))
#Normalize x:
for (j in 1:dim(x)[2]){
  xjbar<-mean(x[,j])
  sigma_j<-sum((x[,j]-xjbar)^2)/(length(x[,j])-1)
  for (i in 1:dim(x)[1]){
    x[i,j]<-(x[i,j]-xjbar)/sqrt(sigma_j)
  }
}

split.select.list <- split.select(x,y,fraction = fraq)
beta_tmp <- split.select.list$beta
lambda <- split.select.list$lambda
split <- split.select.list$split


beta_split <- carve.linear(x,y,split = split, beta = beta_tmp, lambda = lambda,
                         sigma=sigma_squ, normalize_truncation_limits = FALSE)$beta_split

beta_posi <- carve.linear(x,y,split = split, beta = beta_tmp, lambda = lambda,
                           sigma=sigma_squ, normalize_truncation_limits = FALSE)$beta_posi

#Plot the distribution of beta_split, beta_split ~ N(beta^M, sigma^2_MA(X_{B,MA}^T X_{B,MA})^-1)
qqnorm(beta_split)
qqline(beta_split)

length(beta_split)

dim(x[split,])
dim(x[-split,])
sum(beta_tmp!=0)
