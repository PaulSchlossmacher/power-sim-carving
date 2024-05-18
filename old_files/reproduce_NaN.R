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
library(doParallel)
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
#fraq.vec <- c(0.7,0.8,0.9) #to reproduce the error
#fraq.vec <- c(0.5,0.6,0.7,0.8)

#fraq.vec <- c(0.5)
#fraq.vec <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99)
fraq.vec <- c(0.7)
#fraq.vec <- c(0.7)
fraq.vec.Drysdale<-fraq.vec
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
#More active variables than observations in Group B after the split:
#sel.index <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)#active predictors
sel.index <- c(1,5,10,15,20)
beta <- rep(0, p)
beta[sel.index] <- 1
sparsity <- 5
set.seed(42) 
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta
SNR <- 1.713766 # value created for Toeplitz 0.6
sigma_squ <- 2 #Variance 1 instead of 2 before, to make it easier for Lasso to catch the variables
nsim <- 7
sig.level <- 0.05
new_fraq_threshold<-0

total.time <- 0
start.time <- Sys.time()
#Normalize x before starting, y will also be normalized, but at each iteration, as it is always chosen with new noise
for (j in 1:dim(x)[2]){
  xjbar<-mean(x[,j])
  sigma_j<-sum((x[,j]-xjbar)^2)/(length(x[,j])-1)
  for (i in 1:dim(x)[1]){
    x[i,j]<-(x[i,j]-xjbar)/sqrt(sigma_j)
  }
}
conf_matrix <- function(p.vals, sig.level = 0.05, beta){
  #Calculate false positives, true positives, true negatives, false negatives
  H0T_Rej<-sum(p.vals<=sig.level & beta==0)
  H0F_Rej<-sum(p.vals<=sig.level & beta==1)
  H0T_N_Rej<-sum(p.vals>sig.level & beta==0)
  H0F_N_Rej<-sum(p.vals>sig.level & beta==1)
  return (c(H0T_Rej, H0F_Rej, H0T_N_Rej, H0F_N_Rej))
}

plot_conf_matrix <- function(test_res_avg, name, nsim){
  #Create average confusion matrices obtained from one fraction and output to console
  testing_table <- matrix(
    test_res_avg,
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Rejected", "Not Rejected"), c("H0 True", "H0 False"))
  )
  cat("The average confusion matrix of", name ,"p-values over ", nsim, "calculations is: \n")
  print(testing_table)
}


#Initialization of metrics we want to keep track of
f <- length(fraq.vec)
full_test_res_D <- matrix(rep(0,4*f), nrow = f)
full_test_res_C <- matrix(rep(0,4*f), nrow = f)
full_test_res_split <- matrix(rep(0,4*f), nrow = f)
full_test_res_posi <- matrix(rep(0,4*f), nrow = f)
full_power_avg_D <- rep(0,f)
full_power_avg_C <- rep(0,f)
full_power_avg_split <- rep(0,f)
full_power_avg_posi <- rep(0,f)
full_type1_error_avg_D <- rep(0,f)
full_type1_error_avg_C <- rep(0,f)
full_type1_error_avg_split <- rep(0,f)
full_type1_error_avg_posi <- rep(0,f)
full_FWER_D <- rep(0,f)
full_FWER_C <- rep(0,f)
full_FWER_split <- rep(0,f)
full_FWER_posi <- rep(0,f)
#Should count fails of drysdales estimator for a given fraction over nsim rounds
drysdale.fails <- rep(0,f)
RNGkind("L'Ecuyer-CMRG")
#set.seed(42)
state_vector <- c(10407, 291223161, -710388133, 1405033250, -510856512, 1188280219, 202332716)
# Set the RNG state
.Random.seed <- as.integer(state_vector)
print(.Random.seed)
fraq_ind = 1
#for (i in 1:nsim) {
empty_model_C <- FALSE
empty_model_D <- FALSE
counter_new_split <- 0
y <- y.true + sqrt(sigma_squ) * rnorm(n)

#Normalize y:
y<-(y-mean(y))
split.select.list <- split.select(x,y,fraction = fraq.vec[fraq_ind])
beta_tmp <- split.select.list$beta
if(sum(beta_tmp!=0)==0){
  empty_model_C <- TRUE
  empty_model_D <- TRUE
  p_vals_D_fwer <- rep(1,p)
  p_vals_C_fwer <- rep(1,p)
  p_vals_split_fwer <- rep(1,p)
  p_vals_posi_fwer <- rep(1,p)
  print("0 variables where chosen by the lasso, but thats not a problem.")
}
lambda <- split.select.list$lambda
split <- split.select.list$split

# --------- Extra variable selection for beta^Drysdale ------------------

split.select.list_D<-split.select.list
beta_tmp_D<-split.select.list_D$beta
lambda_D <- split.select.list_D$lambda
split_D <- split.select.list_D$split

# While the inverse can not be calculated: 

while(sum(beta_tmp_D!=0)>min(n*fraq.vec.Drysdale[fraq_ind], n*(1-fraq.vec.Drysdale[fraq_ind]))){#CHANGED TO fraq.vec.Drysdale
  #Try new split with less observations for screening:
  
  if (counter_new_split>=new_fraq_threshold){
    fraq.vec.Drysdale[fraq_ind]<-fraq.vec.Drysdale[fraq_ind]-0.025
  }
  
  #CHANGED ALL OF split.selectlist to split.select.list_D
  split.select.list_D <- split.select(x,y,fraction = fraq.vec.Drysdale[fraq_ind])
  beta_tmp_D <- split.select.list_D$beta
  lambda_D <- split.select.list_D$lambda
  split_D <- split.select.list_D$split
  
  counter_new_split<-counter_new_split+1
}

#Note: Since this only applies if the empty model gets selected and the Drysdale controls etc.
#should not kick in here, we don't have to differentiate between different splits here
if(sum(beta_tmp_D!=0)==0){#CHANGED TO beta_tmp_D from beta_tmp, PLEASE CHECK IF THIS IS RIGHT
  empty_model_D <- TRUE
  p_vals_D_fwer <- rep(1,p)
  #p_vals_C_fwer <- rep(1,p) #THIS LINE SHOULD NOT BE NECESSARILY CALLED IF ONLY DRYSDALES MODEL IS EMPTY
  p_vals_split_fwer <- rep(1,p)
  p_vals_posi_fwer <- rep(1,p)
  print("0 variables where chosen by the lasso, but thats not a problem.")
}
#Compute pure p-values from Drysdale's and Christoph's approach
if(!empty_model_C){
  carve_C <- carve.lasso(X = x, y = y, ind = split, beta = beta_tmp, tol.beta = 0, sigma = sigma_squ,
                         lambda = lambda,FWER = FALSE, intercept = FALSE,selected=TRUE, verbose = FALSE)
  p_vals_C_nofwer<-carve_C$pv
  
  #carve_C only returns the p-values of the coefficients determined by the selection event, hence we assign them at the appropriate positions
  p_vals_comp_C<-rep(1,p)
  chosen_C <- which(abs(beta_tmp)>0)
  p_vals_comp_C[chosen_C] <- p_vals_C_nofwer
  
  #Add FWER control with Bonferroni correction
  model.size_C <- length(chosen_C)
  p_vals_C_fwer <- pmin(p_vals_comp_C * model.size_C, 1)
  
  
}
# I decided to include all of posi, split and Drysdale in the case of whether
# Or not beta^Drysdale can be computed, since we work with fraq>0.7 anyways, so
# beta^Posi not working shouldn't really be an issue

if (!empty_model_D){
  carve_D <-carve.linear(x,y,split = split_D, beta = beta_tmp_D, lambda = lambda_D, sigma=sigma_squ)
  p_vals_D_nofwer <- carve_D$pvals
  p_vals_split_nofwer <- beta.split(x, y, split=split_D, beta = beta_tmp_D, sigma=sigma_squ)$pvals_split
  p_vals_posi_nofwer <- beta.posi(x, y, split=split_D, beta = beta_tmp_D,lambda=lambda_D, sigma=sigma_squ)$pvals#CHANGED LAMBDA TO LAMBDA_D
  
  chosen_D <- which(abs(beta_tmp_D)>0)
  model.size_D<- length(chosen_D)
  p_vals_D_fwer <- pmin(p_vals_D_nofwer * model.size_D, 1)
  p_vals_split_fwer <- pmin(p_vals_split_nofwer*model.size_D,1)
  p_vals_posi_fwer <- pmin(p_vals_posi_nofwer*model.size_D,1)
  
}
  # if is.na(p_vals_posi_fwer){
  #   browser()
  # }
#}
