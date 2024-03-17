#rm(list=ls())#leave this if we want to empty all global variables

#TODO:add these repositories to github and adapt the paths to work for all of us without modification
hdi_adjustments_path <- "C:/Users/Filip/Documents/ETH/Master/Semester 2/Semesterarbeit/Multicarving Christoph/inference/hdi_adjustments.R"
carving_path <- "C:/Users/Filip/Documents/ETH/Master/Semester 2/Semesterarbeit/Multicarving Christoph/inference/carving.R"
sample_from_truncated_path <- "C:/Users/Filip/Documents/ETH/Master/Semester 2/Semesterarbeit/Multicarving Christoph/inference/sample_from_truncated.R"
tryCatchWE_path <- "C:/Users/Filip/Documents/ETH/Master/Semester 2/Semesterarbeit/Multicarving Christoph/inference/tryCatch-W-E.R"
SNTN_distribution <- "C:/Users/Filip/Documents/ETH/Master/Semester 2/Semesterarbeit/Semester-Project-Multicarving/SNTN_distribution.R"
split_select_function <- "C:/Users/Filip/Documents/ETH/Master/Semester 2/Semesterarbeit/Semester-Project-Multicarving/split_select.R"


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
source(SNTN_distribution)
source(split_select_function)


# toeplitz
n <- 100
p <- 200
rho <- 0.6
fraq = 0.9
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
sel.index <- c(1, 5, 10, 15, 20)#active predictors
ind <- sel.index
beta <- rep(0, p)#initialize beta as all zeros
beta[sel.index] <- 1#put ones at active predictor positions
sparsity <- 5
set.seed(42) # to make different methods comparable, fix the x-matrix
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta
SNR <- 1.713766 # value created for Toeplitz 0.6
sigma <- 2
y <- y.true + sigma * rnorm(n)

pvals.v <- rep(1,p)
sel.models <- logical(p)#logical vector initialized with False, this will indicate which predictors are chosen in the model


#make a split that yields valid constraints and perform model selection on it
split_select_list <- split_select(x,y,fraction = fraq)
beta <- split_select_list$beta
split <- split_select_list$split
x_a <- x[split, ]
y_a <- y[split,]
x_b <- x[-split, ]
y_b <- y[-split]

chosen <-  which(abs(beta) > 0) # selected variables
s <- length(chosen)
if (s == 0) {
  return(NULL)
}
#Did not account for intercept until now, need to find out how to control selection with or without intercept
# if (intercept && s == 1) {
#   return(NULL)
# }

#extract active variables from both splits
x_Ma <- x_a[, chosen]
x_Mb <- x_b[, chosen]

#computes the moore penrose inverse of active variables in both splits
x_Ma.i <- ginv(x_Ma)
x_Mb.i <- ginv(x_Mb)

beta_carve <- 1/n*(fraq*x_Ma.i%*%y_a + (1-fraq)*x_Mb.i%*%y_b)

#inverse of information matrix in linear regression: (X^T*X)^-1
inv.info_a <- tcrossprod(x_Ma.i) #tcrossprod(X,X) does X%*%X^T slightly faster than using the transpose explicitly
inv.info_b <- tcrossprod(x_Mb.i)

#should maybe use the tools from this page somehow to obtain v_lo and v_up, 
#otherwise we cannot determine omega and delta from lemma 3.2 in drysdale:
#https://cran.r-project.org/web/packages/selectiveInference/selectiveInference.pdf
#i did not yet realize how to obtain these truncation limits from christophs code

#DUMMY TEST for sntn_cdf, exchange these values with the corresponding ones for the distribution of beta_carve
# mu1 <- 0
# tau1 <- 1
# mu2 <- 1
# tau2 <- 2
# a <- 2
# b <- 5
# c1 <- 0.5
# c2 <- 0.5
# # Test value for z
# z <- 1.5
# #call sntn
# sntn_cdf <- SNTN_CDF(z, mu1, tau1, mu2, tau2, a, b, c1, c2)
# # Print the result5
# print(sntn_cdf)
