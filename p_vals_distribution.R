#This file creates a visualization of the distribution of p-values of truly inactive coefficients in beta_comb 
#under screening and without screening
#Clear all variables
rm(list = ls())


#Local & user specific path
Local_path<-getwd()
hdi_adjustments_path<-paste(Local_path, "/multicarving_paper/inference/hdi_adjustments.R", sep="")
carving_path<-paste(Local_path, "/multicarving_paper/inference/carving.R", sep="")
sample_from_truncated_path<-paste(Local_path, "/multicarving_paper/inference/sample_from_truncated.R", sep="")
tryCatchWE_path<-paste(Local_path, "/multicarving_paper/inference/tryCatch-W-E.R", sep="")

#Different paths here, because they're in a different directory
SNTN_distribution_path<-paste(Local_path, "/SNTN_distribution.R", sep="")
split_select_function_path<-paste(Local_path, "/split_select.R", sep="")
carve_combined_path<-paste(Local_path, "/carve_combined.R", sep="")


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
library(ggplot2)
library(dplyr)
library(gridExtra)

source(hdi_adjustments_path)
source(carving_path)
source(sample_from_truncated_path)
source(tryCatchWE_path)
source(SNTN_distribution_path)
source(split_select_function_path)
source(carve_combined_path)



#-------------------- Toeplitz Carving simulation from Christoph ----------------------
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
SNR = 2
sigma_squ = drop(var(y.true)) / SNR

#Normalize x before starting, y will also be normalized, but at each iteration, as it is always chosen with new noise
for (j in 1:dim(x)[2]){
  xjbar<-mean(x[,j])
  sigma_j<-sum((x[,j]-xjbar)^2)/(length(x[,j])-1)
  for (i in 1:dim(x)[1]){
    x[i,j]<-(x[i,j]-xjbar)/sqrt(sigma_j)
  }
}

counter <- 0
screening <- c()

#This will contain all predictor p-values, which are truly inactive but where chosen in the selection event
#Just collect all of them in one vector, as we dont care in which simulation round we obtained them
p_vals_screen <- c()
p_vals_noscreen <- c()
target_number <- 10000
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
      stop("Tried to many selection events and not one of them was conformable for carve.comb")
    }
    select.again <- FALSE
    set.seed(counter)
    counter <- counter + 1
    y <- y.true + sqrt(sigma_squ) * rnorm(n)
    y <- y - mean(y)
    split.select.list <- split.select(x,y,fraction = fraq)
    beta_tmp <- split.select.list$beta
    
    if(sum(beta_tmp!=0)==0){
      select.again <- TRUE
      print("0 variables where chosen by the lasso, repeating selection")
    }
    lambda <- split.select.list$lambda
    split <- split.select.list$split
    if(sum(beta_tmp!=0)>min(n*fraq, n*(1-fraq))){
      select.again <- TRUE
      select.again.counter <- select.again.counter + 1
      print("Need to split again because we selected more variables than carve.comb can handle")
    }
    
  }
  p_vals_D <- carve.comb(x,y,split = split, beta = beta_tmp, lambda = lambda,sigma_squ=sigma_squ)$pvals
  
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


print(length(p_vals_screen))
print(length(p_vals_noscreen))

p_vals_screen_df <- data.frame(p_values = p_vals_screen, Type = "Screening")
p_vals_noscreen_df <- data.frame(p_values = p_vals_noscreen, Type = "No Screening")
p_vals_df <- rbind(p_vals_screen_df, p_vals_noscreen_df)

plot_screen <- ggplot(p_vals_screen_df, aes(x = p_values)) +
  stat_ecdf(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = expression("Screening"),
    x = "x",
    y = expression(F[n](x))
  ) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(0,1))+
  scale_y_continuous(limits=c(0,1))

plot_noscreen <- ggplot(p_vals_noscreen_df, aes(x = p_values)) +
  stat_ecdf(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = expression("No screening"),
    x = "x",
    y = expression(F[n](x))
  ) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(0,1))+
  scale_y_continuous(limits=c(0,1))


#print the obtained plots side by side
grid.arrange(plot_screen, plot_noscreen, ncol = 2)
#Safe them as a single image
g <- arrangeGrob(plot_screen, plot_noscreen, ncol = 2)
ggsave("ecdf_plots_screen_and_noscreen.png", g, width = 10, height = 5)

#Optional: plot a histogram of p-values under screening
#Optional: plot a histogram of p-values under screening
# hist(p_vals_screen_df$p_values,
#      breaks = seq(0, 1, by = 0.1),
#      main = expression("Histogram of" ~ italic("p") ~ "-values under Screening"),
#      xlab = "p-values",
#      ylab = "Frequency",
#      col = "darkblue",
#      border = "black",
#      xlim = c(0, 1),
#      freq = TRUE)


