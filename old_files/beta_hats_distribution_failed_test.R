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
n <- 100
p <- 200
rho <- 0.6
fraq = 0.7
#fraq.vec <- c(0.7)
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
#More active variables than observations in Group B after the split:
act.index <- c(1, 5, 10, 15, 20)#active predictors
beta <- rep(0, p)
beta[act.index] <- 1
set.seed(42)
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta
sigma <- 1.2

# ------------- Different set of parameters for testing the behaviour: -------------
# 
# n <- 100
# p <- 200
# rho <- 0.7
# fraq = 0.7
# #fraq.vec <- c(0.7)
# #toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
# Cov <- toeplitz(rho ^ (seq(0, p - 1)))
# #More active variables than observations in Group B after the split:
# act.index <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)#active predictors
# beta <- rep(0, p)
# beta[act.index] <- 1
# #sparsity <- 5
# set.seed(1)
# x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
# y.true <- x %*% beta
# SNR <- 1.713766 # value created for Toeplitz 0.6
# sigma <- 2

# --------------------------------------------------------

#Normalize x before starting, y will also be normalized, but at each iteration, as it is always chosen with new noise
for (j in 1:dim(x)[2]){
  xjbar<-mean(x[,j])
  sigma_j<-sum((x[,j]-xjbar)^2)/(length(x[,j])-1)
  for (i in 1:dim(x)[1]){
    x[i,j]<-(x[i,j]-xjbar)/sqrt(sigma_j)
  }
}

#nsim = 20
counter <- 1
screening <- c()

#This will contain all beta_split and beta_posi coefficients
#Just collect all of them in one vector, as we dont care in which simulation round we obtained them
beta_split_list <- c()
beta_posi_list <- c()
target_number <- 3000
# We could probably choose a smaller size as well, but this shouldn't be too
# computationally expensive anyways
# Note: set.seed(42) from line 57 for setting these seeds
possible_seeds<-sample(x=1:100000, replace=FALSE, size=target_number)

beta_split_count <-  0
beta_posi_count <- 0
rounds <- 0
#Note: beta_split and beta_posi should have the same length anyways, no need for seperate counters
while (beta_split_count < target_number){
  rounds <- rounds + 1
  #get different selection events
  select.again <- TRUE
  select.again.counter = 0
  while(select.again){
    if (select.again.counter > 50){
      stop("Tried to many selection events and not one of them was conformable for beta_Drysdale")
    }
    select.again <- FALSE
    set.seed(possible_seeds[counter])
    counter <- counter + 1
    y <- y.true + sigma * rnorm(n)
    #Normalize y:
    y<-(y-mean(y))
    
    split.select.list <- split.select(x,y,fraction = fraq)
    beta_tmp <- split.select.list$beta
    chosen_indices<-which(beta_tmp!=0)
    
    if(sum(beta_tmp!=0)==0){
      select.again <- TRUE
      print("0 variables where chosen by the lasso, repeating selection")
    }

    lambda <- split.select.list$lambda
    split <- split.select.list$split
    if(sum(beta_tmp!=0)>min(n*fraq, n*(1-fraq))){
      select.again <- TRUE
      select.again.counter <- select.again.counter + 1
      print("Need to split again because we selected more variables than beta_D can handle")
    }
    
  }
  
  beta_split<-carve.linear(x,y,split = split, beta = beta_tmp, lambda = lambda,
                           sigma=sigma, normalize_truncation_limits = FALSE)$beta_split
  beta_posi<-carve.linear(x,y,split = split, beta = beta_tmp, lambda = lambda,
                           sigma=sigma, normalize_truncation_limits = FALSE)$beta_posi
  
  
  sel.index <- which(beta_tmp != 0)
  #check screening condition
  if (all(act.index %in% sel.index)){
    beta_split_list<-c(beta_split_list, beta_split)
    beta_posi_list<-c(beta_posi_list, beta_posi)    
  
    beta_split_count<-length(beta_split_list)
    beta_posi_count<-length(beta_posi_list)

    screening <- c(screening, TRUE)
    
  }
  #if screening not successful skip the round
  else{
    screening <- c(screening, FALSE)
    next
  }

  cat("Current length of beta_hat values is:", beta_split_count, "\n")

  
}

qqnorm(beta_split_list)
qqline(beta_split_list)

plot(density(beta_split_list), main="Empirical PDF of beta_split", xlab="Values", ylab="Density")

# After getting the results above, I wanted to try the same for values, for which I know
# that they are normally distributed. However it turns out that just concatenating 
# normally distributed vectors does not lead to a normal distribution (could have thought
# of that sooner tbh)

test<-c(rnorm(5000, mean=0, sd=1), rnorm(5000, mean=4, sd=5), rnorm(5000, mean=-4, sd=8))

plot(density(test), main="Empirical PDF Test", xlab="Values", ylab="Density")

# Therefore one should disregard the results of this file. However I did not want to delete it either, just to be safe
