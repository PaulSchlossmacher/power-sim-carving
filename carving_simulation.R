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


source(hdi_adjustments_path)
source(carving_path)
source(sample_from_truncated_path)
source(tryCatchWE_path)
source(SNTN_distribution_path)
source(split_select_function_path)
source(carve_linear_path)

#-------------------- Toeplitz Carving simulation from Christoph ----------------------
n <- 100
p <- 200
rho <- 0.6
fraq = 0.9
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
sel.index <- c(1, 5, 10, 15, 20)#active predictors
beta <- rep(0, p)#initialize beta as all zeros
beta[sel.index] <- 1#put ones at active predictor positions
sparsity <- 5
set.seed(42) # to make different methods comparable, fix the x-matrix
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta
SNR <- 1.713766 # value created for Toeplitz 0.6
sigma <- 2
y <- y.true + sigma * rnorm(n)


#Here we run the simulation with beta_Carve^Drysdale
carve_D <-carve.linear(x,y,fraq,sigma=sigma)
split <- carve_D$split
beta_tmp <- carve_D$beta
lambda <- carve_D$lambda

#We get some warnings for hamiltonian sampler
carve_C <- carve.lasso(X = x, y = y, ind = split, beta = beta_tmp, tol.beta = 0, sigma = sigma,
                             lambda = lambda, intercept = FALSE,selected=TRUE, verbose = TRUE)