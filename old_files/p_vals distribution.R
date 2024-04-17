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
fraq = 0.7
#fraq.vec <- c(0.7)
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
sigma <- 1.9 #Variance 1.9 instead of 2 before, to make screening about as often successful as not

#nsim = 20
counter <- 0
screening <- c()

#This will contain all predictor p-values, which are truly inactive but where chosen in the selection event
#Just collect all of them in one vector, as we dont care in which simulation round we obtained them
p_vals_screen <- c()
p_vals_noscreen <- c()
target_number <- 200
p_val_screen_count <-  0
p_val_noscreen_count <- 0
rounds <- 0
while (p_val_screen_count < target_number || p_val_noscreen_count < target_number){
  rounds <- rounds + 1
  #get different selection events
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
  
  p_vals_D <- carve.linear(x,y,split = split, beta = beta_tmp, lambda = lambda, fraction = fraq,sigma=sigma)
  
  sel.index <- which(beta_tmp != 0)
  #check screening condition
  if (all(act.index %in% sel.index)){
    if(p_val_screen_count > target_number){
      next
    }
    screening <- c(screening, TRUE)
    p_vals_screen <- c(p_vals_screen, p_vals_D[sel.index][!(sel.index %in% act.index)])
    p_val_screen_count <- length(p_vals_screen)
  }
  else{
    if(p_val_noscreen_count > target_number){
      next
    }
    p_vals_noscreen <- c(p_vals_noscreen, p_vals_D[sel.index][!(sel.index %in% act.index)])
    screening <- c(screening, FALSE)
    p_val_noscreen_count <- length(p_vals_noscreen)
  }
  print(min(p_val_screen_count,p_val_noscreen_count))
  
 }
cat("We had ", sum(screening), " successful screenings out of ", rounds, " simulations.")

# par(mfrow = c(1, 2))
# #Create empty plot to visualize 1dim distribution of p-values
# plot(x = NULL, y = NULL, xlim = c(0, 1), ylim = c(0, 1), 
#      main = "p-values under screening", xlab = "Values", ylab = "")
# 
# # Add points to the plot at a fixed y-coordinate
# points(p_vals_screen, rep(0.5, length(p_vals_screen)), pch = 16, col = "skyblue", cex = 0.7)
# 
# plot(x = NULL, y = NULL, xlim = c(0, 1), ylim = c(0, 1), 
#      main = "p-values without screening", xlab = "Values", ylab = "")
# 
# # Add points to the plot at a fixed y-coordinate
# points(p_vals_noscreen, rep(0.5, length(p_vals_noscreen)), pch = 16, col = "skyblue", cex = 0.7)

print(length(p_vals_screen))
print(length(p_vals_noscreen))

#Create and save QQ-plots to see how well the p-values match a uniform distribution
quantiles_uniform <- runif(target_number)
png("QQ_plots.png",width = 800, height = 400)
par(mfrow = c(1, 2))
qqplot(quantiles_uniform, p_vals_screen[1:target_number],
       xlab = "Theoretical Quantiles", ylab = "p_vals_screen", main = "QQ-Plot: P-values under screening")
abline(0, 1, col = "red")


qqplot(quantiles_uniform, p_vals_noscreen[1:target_number],
       xlab = "Theoretical Quantiles", ylab = "p_vals_noscreen", main = "QQ-Plot: P-values without screening")
abline(0, 1, col = "red")
dev.off()

#Plot a histogram of p_values under screening
png("histogram.png",width = 350, height = 300)
hist(p_vals_screen, main = "Histogram of p-values under screening")
dev.off()
