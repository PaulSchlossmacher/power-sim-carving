#Trying to create a more extreme Toeplitz example with less noise to 
#encourage Lasso to use more variables:
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
fraq.vec <- c(0.7)
#fraq.vec <- c(0.5,0.6,0.7)
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
sigma_squ <- 1 #Variance 1 instead of 2 before, to make it easier for Lasso to catch the variables

nsim <- 3
f <- length(fraq.vec)
full_test_res_D <- matrix(rep(0,4*f), nrow = f)
full_test_res_C <- matrix(rep(0,4*f), nrow = f)
full_power_avg_D <- rep(0,f)
full_power_avg_C <- rep(0,f)
full_type1_error_avg_D <- rep(0,f)
full_type1_error_avg_C <- rep(0,f)

#Should count fails of drysdales estimator for a given fraction over nsim rounds
drysdale.fails <- rep(0,f)

for (fraq_ind in 1:f){


  test_res_D <- rep(0,4)
  test_res_C <- rep(0,4)
  counter <- 0
  powers_D <- rep(0,nsim)
  powers_C <- rep(0,nsim)
  type1_error_D <- rep(0,nsim)
  type1_error_C <- rep(0,nsim)
  
  
  for (i in 1:nsim){
    print(i)
    #get different selection events
    select.again <- TRUE
    select.again.counter = 0
    while(select.again){
      if (select.again.counter > 200){
        stop("Tried to many selection events and not one of them was conformable for beta_Drysdale")
      }
      select.again <- FALSE
      set.seed(counter)
      counter <- counter + 1
      y <- y.true + sqrt(sigma_squ) * rnorm(n)
      split.select.list <- split.select(x,y,fraction = fraq.vec[fraq_ind])
      beta_tmp <- split.select.list$beta
      if(sum(beta_tmp!=0)==0){
        #TODO:set entries of power and type 1 error vectors to NA if no variables where selected?? and skip round with "next"
        select.again <- TRUE
        print("0 variables where chosen by the lasso, repeating selection")
      }
      lambda <- split.select.list$lambda
      split <- split.select.list$split
      if(sum(beta_tmp!=0)>min(n*fraq.vec[fraq_ind], n*(1-fraq.vec[fraq_ind]))){
        select.again <- TRUE
        select.again.counter <- select.again.counter + 1
        print("Need to split again because we selected more variables than beta_D can handle")
      }
    }
    
    #print("calculating Drysdales p-values")
    carve_D <-carve.linear(x,y,split = split, beta = beta_tmp, lambda = lambda, sigma=sigma_squ)
    p_vals_D <- carve_D$pvals
    
    #false positives, true positives, true negatives, false negatives for Drysdales p-values
    H0T_Rej_D<-sum(p_vals_D<=0.05 & beta==0)
    H0F_Rej_D<-sum(p_vals_D<=0.05 & beta==1)
    H0T_N_Rej_D<-sum(p_vals_D>0.05 & beta==0)
    H0F_N_Rej_D<-sum(p_vals_D>0.05 & beta==1)
    
    #print("calculating Christophs p-values")
    carve_C <- carve.lasso(X = x, y = y, ind = split, beta = beta_tmp, tol.beta = 0, sigma = sigma_squ,
                           lambda = lambda, intercept = FALSE,selected=TRUE, verbose = FALSE)
    p_vals_C<-carve_C$pv
    p_vals_comp_C<-rep(1,p)
    chosen_C <- which(abs(beta_tmp)>0)
    p_vals_comp_C[chosen_C] <- p_vals_C
    
    #false positives, true positives, true negatives, false negatives for Drysdales p-values
    H0T_Rej_C<-sum(p_vals_comp_C<=0.05 & beta==0)
    H0F_Rej_C<-sum(p_vals_comp_C<=0.05 & beta==1)
    H0T_N_Rej_C<-sum(p_vals_comp_C>0.05 & beta==0)
    H0F_N_Rej_C<-sum(p_vals_comp_C>0.05 & beta==1)

    #Collecting terms
    test_res_D<-test_res_D + c(H0T_Rej_D,H0F_Rej_D,H0T_N_Rej_D,H0F_N_Rej_D)
    test_res_C<- test_res_C + c(H0T_Rej_C,H0F_Rej_C,H0T_N_Rej_C,H0F_N_Rej_C)
    #calculating power and type1 error for 1 round in simulation
    powers_D[i] <- H0F_Rej_D/(length(sel.index))
    powers_C[i] <- H0F_Rej_C/(length(sel.index))
    type1_error_D[i] <- H0T_Rej_D/(p-length(sel.index))
    type1_error_C[i] <- H0T_Rej_C/(p-length(sel.index))
    
  }
  
  #calculating averages over all rounds of one fraction
  test_res_D_avg <- test_res_D/nsim
  test_res_C_avg <- test_res_C/nsim
  power_avg_D <- mean(powers_D)
  power_avg_C <- mean(powers_C)
  type1_error_avg_D <- mean(type1_error_D)
  type1_error_avg_C <- mean(type1_error_C)
  
  #Creating average confusion matrices obtained from one fraction
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
  #Printing results of one fraction to console
  cat("Results for fraction", fraq.vec[fraq_ind], ":")
  cat("The average confusion matrix of Drydales p-values over ", nsim, "calculations is: \n" )
  print(Testing_table_D_avg)
  cat("The average confusion matrix of Christoph's p-values over ", nsim, "calculations is: \n" )
  print(Testing_table_C_avg)
  
  cat("The average power of Drysdales p-values:", power_avg_D)
  cat("The average power of Christophs p-values:", power_avg_C)
  cat("The average type 1 error of Drysdales p-values:", type1_error_avg_D)
  cat("The average type 1 error of Christophs p-values:", type1_error_avg_C)
  
  
  #Store everything to create plots later
  full_test_res_D[fraq_ind, ] <- test_res_D_avg
  full_test_res_C[fraq_ind, ] <- test_res_C_avg
  full_power_avg_D[fraq_ind] <- power_avg_D
  full_power_avg_C[fraq_ind] <- power_avg_C
  full_type1_error_avg_D[fraq_ind] <- type1_error_avg_D
  full_type1_error_avg_C[fraq_ind] <- type1_error_avg_C
  #collecting fails of drysdales estimator of one fraction 
  drysdale.fails[fraq_ind] <- counter - nsim
}

#save.image(file='myEnvironment.RData')
#load('myEnvironment.RData')



# --------------- Create plots --------------

data_Power <- data.frame(
  Fraq=fraq.vec,
  "Avg Power Christoph" = full_power_avg_C,
  "Avg Power Drysdale" = full_power_avg_D
)

data_Power_long <- tidyr::gather(data_Power, "Type", "Value", -Fraq)

PowerPlot<-ggplot(data_Power_long, aes(x = Fraq, y = Value, color = Type)) +
  geom_line() +
  labs(title = "Average Power",
       x = "Fraq", y = "Value") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5))

ggsave("PowerPlot.png", plot = PowerPlot, width = 8, height = 6,
       units = "in", dpi = 300, bg = "#F0F0F0")


data_TypeI <- data.frame(
  Fraq=fraq.vec,
  "Avg Type I Error rate Christoph" = full_type1_error_avg_C,
  "Avg Type I Error rate Drysdale" = full_type1_error_avg_D
)

data_TypeI_long <- tidyr::gather(data_TypeI, "Type", "Value", -Fraq)

TypeIPlot<-ggplot(data_TypeI_long, aes(x = Fraq, y = Value, color = Type)) +
  geom_line() +
  labs(title = "Average Type I Error Rate",
       x = "Fraq", y = "Value") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5))

ggsave("TypeIPlot.png", plot = TypeIPlot, width = 8, height = 6,
       units = "in", dpi = 300, bg = "#F0F0F0")
