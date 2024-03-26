#Trying to create a more extreme Toeplitz example with less noise to 
#encourage Lasso to use more variables:



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
fraq = 0.6
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
#More active variables than observations in Group B after the split:
sel.index <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)#active predictors
beta <- rep(0, p)
beta[sel.index] <- 1
sparsity <- 5
set.seed(42) 
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta
SNR <- 1.713766 # value created for Toeplitz 0.6
sigma <- 1 #Variance 1 instead of 2 before, to make it easier for Lasso to catch the variables
y <- y.true + sigma * rnorm(n)

# Tried normal glmnet Lasso just to see whether screening works on the full sample:
# It does.
# cv_model <- cv.glmnet(x, y, alpha = 1)
# plot(cv_model)
# coef(cv_model)


#Here we run the simulation with beta_Carve^Drysdale
#carve_D <-carve.linear(x,y,fraq,sigma=sigma)

#We get an error 

#Not necessary to run those, if we get an error before:
# split <- carve_D$split
# beta_tmp <- carve_D$beta
# lambda <- carve_D$lambda

#We get an error message

split.select.list <- split.select(x,y,fraction = fraq)
beta_tmp <- split.select.list$beta
lambda <- split.select.list$lambda
split <- split.select.list$split


carve_C <- carve.lasso(X = x, y = y, ind = split, beta = beta_tmp, tol.beta = 0, sigma = sigma,
                       lambda = lambda, intercept = FALSE,selected=TRUE, verbose = TRUE)
carve_C
p_vals_C<-carve_C$pv
#For carve_C we still get p-values


#Let's try to make beta^Drysdale work, by choosing a different split:
sum(beta_tmp!=0)

#We have 37 non-0 beta coefficients. So let's go for a 60-40 split instead of 90-10 before:

carve_D <-carve.linear(x,y,fraction=0.6,sigma=sigma)
beta_tmp2 <- carve_D$beta
p_vals_D <- carve_D$pvals
p_vals_comp_D <- rep(1,p)
chosen_D <- which(abs(beta_tmp2)>0)
p_vals_comp_D[chosen_D] <- p_vals_D


#Now we don't get an error anymore and we can calculate the p-values:


#Create a sort of confusion matrix for carve_C:
#True Positives:
#We first have to create a new vector of p-values in R^p, i.e. R^200, where all
#variables that weren't selected get assigned a p-value of 1:
#Filip: maybe better readability without for loop
p_vals_comp_C<-rep(1,p)
chosen_C <- which(abs(beta_tmp)>0)
p_vals_comp_C[chosen_C] <- p_vals_C
# Counter=1
# for (i in 1:p){
#   if (beta_tmp[i]!=0){
#     p_vals_comp_C[i]<-p_vals_C[Counter]
#     Counter=Counter+1
#   }
# }
#Upon first viewing this seems to work

#Question to Paul: What happens here?
#prd<-ifelse(blabla>05., "Yes", "No")
#table(prd, d$response)

#Do it with level of significance 0.05, but maybe sth else is smarter?
#Multiple testing issue?

H0T_Rej_C<-sum(p_vals_comp_C<=0.05 & beta==0)
H0F_Rej_C<-sum(p_vals_comp_C<=0.05 & beta==1)
H0T_N_Rej_C<-sum(p_vals_comp_C>0.05 & beta==0)
H0F_N_Rej_C<-sum(p_vals_comp_C>0.05 & beta==1)



Testing_table_C <- matrix(
  c(H0T_Rej_C, H0F_Rej_C, H0T_N_Rej_C, H0F_N_Rej_C),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("Rejected", "Not Rejected"), c("H0 True", "H0 False"))
)

print(Testing_table_C)

# p_vals_comp_C[65] with a p-val of 0.097 is the only one we don't catch


#Those are Christophs results, what's left is to do the same for Drysdale.
#Also this is of course a very small sample size for comparison and "power studies"
#Following Pauls procedure, here is the same for Drysdales beta_carve
H0T_Rej_D<-sum(p_vals_comp_D<=0.05 & beta==0)
H0F_Rej_D<-sum(p_vals_comp_D<=0.05 & beta==1)
H0T_N_Rej_D<-sum(p_vals_comp_D>0.05 & beta==0)
H0F_N_Rej_D<-sum(p_vals_comp_D>0.05 & beta==1)



Testing_table_D <- matrix(
  c(H0T_Rej_D, H0F_Rej_D, H0T_N_Rej_D, H0F_N_Rej_D),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("Rejected", "Not Rejected"), c("H0 True", "H0 False"))
)

print(Testing_table_D)

#Maybe we should try running a simulation that does all of the above for example a 100 times?