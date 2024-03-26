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
fraq = 0.7
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


# Tried normal glmnet Lasso just to see whether screening works on the full sample:
# It does.
# cv_model <- cv.glmnet(x, y, alpha = 1)
# plot(cv_model)
# coef(cv_model)

test_res_D <- rep(0,4)
test_res_C <- rep(0,4)
counter <- 1
nsim <- 3
powers_D <- rep(0,nsim)
powers_C <- rep(0,nsim)
type1_error_D <- rep(0,nsim)
type1_error_C <- rep(0,nsim)


for (i in 1:nsim){
  #get different selection events
  select.again <- TRUE
  while(select.again){
    if (counter > 50){
      stop("Tried to many selection events and not one of them was conformable for beta_Drysdale")
    }
    select.again <- FALSE
    set.seed(counter)
    counter <- counter + 1
    y <- y.true + sigma * rnorm(n)
    split.select.list <- split.select(x,y,fraction = fraq)
    beta_tmp <- split.select.list$beta
    lambda <- split.select.list$lambda
    split <- split.select.list$split
    if(sum(beta_tmp!=0)>min(n*fraq, n*(1-fraq))){
      select.again <- TRUE
      print("Need to split again because we selected more variables than beta_D can handle")
    }
  }
  
  print("calculating Drysdales p-values")
  p_vals_D <-carve.linear(x,y,split = split, beta = beta_tmp, lambda = lambda, fraction = fraq,sigma=sigma)
  
  print("calculating Christophs p-values")
  carve_C <- carve.lasso(X = x, y = y, ind = split, beta = beta_tmp, tol.beta = 0, sigma = sigma,
                         lambda = lambda, intercept = FALSE,selected=TRUE, verbose = FALSE)
  p_vals_C<-carve_C$pv
  p_vals_comp_C<-rep(1,p)
  chosen_C <- which(abs(beta_tmp)>0)
  p_vals_comp_C[chosen_C] <- p_vals_C
  
  
  
  #Create a sort of confusion matrix for carve_C:
  #True Positives:
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
  
  
  
  # p_vals_comp_C[65] with a p-val of 0.097 is the only one we don't catch
  
  
  #Those are Christophs results, what's left is to do the same for Drysdale.
  #Also this is of course a very small sample size for comparison and "power studies"
  #Following Pauls procedure, here is the same for Drysdales beta_carve
  H0T_Rej_D<-sum(p_vals_D<=0.05 & beta==0)
  H0F_Rej_D<-sum(p_vals_D<=0.05 & beta==1)
  H0T_N_Rej_D<-sum(p_vals_D>0.05 & beta==0)
  H0F_N_Rej_D<-sum(p_vals_D>0.05 & beta==1)
  
  
  
  Testing_table_D <- matrix(
    c(H0T_Rej_D, H0F_Rej_D, H0T_N_Rej_D, H0F_N_Rej_D),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Rejected", "Not Rejected"), c("H0 True", "H0 False"))
  )
  test_res_D<-test_res_D + c(H0T_Rej_D,H0F_Rej_D,H0T_N_Rej_D,H0F_N_Rej_D)
  test_res_C<- test_res_C + c(H0T_Rej_C,H0F_Rej_C,H0T_N_Rej_C,H0F_N_Rej_C)
  
  powers_D[i] <- H0F_Rej_D/(length(sel.index))
  powers_C[i] <- H0F_Rej_C/(length(sel.index))
  type1_error_D[i] <- H0T_Rej_D/(p-length(sel.index))
  type1_error_C[i] <- H0T_Rej_C/(p-length(sel.index))
  
  
  
  print(Testing_table_D)
  print(Testing_table_C)
}

test_res_D_avg <- test_res_D/nsim
test_res_C_avg <- test_res_C/nsim

power_avg_D <- mean(powers_D)
power_avg_C <- mean(powers_C)
type1_error_avg_D <- mean(type1_error_D)
type1_error_avg_C <- mean(type1_error_C)

Testing_table_D_avg <- matrix(
  test_res_D_avg,
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("Rejected", "Not Rejected"), c("H0 True", "H0 False"))
)

Testing_table_C_avg <- matrix(
  test_res_C_avg,
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("Rejected", "Not Rejected"), c("H0 True", "H0 False"))
)
cat("The average confusion matrix of Drydales p-values over ", nsim, "calculations is: \n" )
print(Testing_table_D_avg)
cat("The average confusion matrix of Christoph's p-values over ", nsim, "calculations is: \n" )
print(Testing_table_C_avg)

cat("The average power of Drysdales p-values:", power_avg_D)
cat("The average power of Christophs p-values:", power_avg_C)
cat("The average type 1 error of Drysdales p-values:", type1_error_avg_D)
cat("The average type 1 error of Christophs p-values:", type1_error_avg_C)


#Maybe we should try running a simulation that does all of the above for example a 100 times?