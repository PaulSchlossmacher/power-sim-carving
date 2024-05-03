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

fraq.vec <- c(0.5, 0.95)
#fraq.vec <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99)

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
nsim <- 10
sig.level <- 0.05
new_fraq_threshold<-0
fraq.vec.Drysdale<-fraq.vec

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

progress <- function(n, tag) {
  mod <- 12
  if (n %% mod == 0 ) {
    cat(sprintf('tasks completed: %d; tag: %d\n', n, tag))
  }
  if (n %% mod == 0 ) {
    toc()
    tic()
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
# To keep track of the average fraq used for beta^Drysdale over all simulation runs:
avg_fraq.vec.Drysdale<-rep(0,f)

RNGkind("L'Ecuyer-CMRG")
set.seed(42)
seed.vec <- sample(1:10000, length(fraq.vec))
seed.n <- 0

for(fraq_ind in  1:f){
  seed.n <- seed.n + 1
  set.seed(seed.vec[seed.n])
  
  opts <- list(progress = progress)
  cl<-makeSOCKcluster(12) 
  rseed <- seed.vec[seed.n]
  clusterSetRNGStream(cl, iseed = rseed) #make things reproducible
  registerDoSNOW(cl)
  tic()
  #start parallel computation
  results <- foreach(i = 1:nsim,.combine = 'rbind', .multicombine = TRUE, 
                     .packages = c("MASS", "mvtnorm", "glmnet", "Matrix", "tictoc", 
                                   "hdi", "selectiveInference", "truncnorm"), .options.snow = opts) %dorng%{
                                     #Initialise fraq.vec.Drysdale
                                     fraq.vec.Drysdale <- fraq.vec
                                     
                                     
                                     #get different selection events
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
                                       print("0 variables where chosen by the lasso, but thats not a problem.t")
                                     }
                                     lambda <- split.select.list$lambda
                                     split <- split.select.list$split
                                     
                                     # --------- Extra variable selection for beta^Drysdale ------------------
                                     
                                     split.select.list_D<-split.select.list
                                     beta_tmp_D<-split.select.list_D$beta
                                     lambda_D <- split.select.list_D$lambda
                                     split_D <- split.select.list_D$split
                                     
                                     # While the inverse can not be calculated: 
                                     
                                     #Note: We discarded the original condition with the minimum, because the Lasso should only select
                                     # s<=min{n_A,p} variables anyways.
                                     #while(sum(beta_tmp_D!=0)>min(n*fraq.vec.Drysdale[fraq_ind], n*(1-fraq.vec.Drysdale[fraq_ind]))){
                                     while(sum(beta_tmp_D!=0)>n*(1-fraq.vec.Drysdale[fraq_ind])){
                                       #Try new split with less observations for screening:
                                       
                                       if (counter_new_split>=new_fraq_threshold){
                                         fraq.vec.Drysdale[fraq_ind]<-fraq.vec.Drysdale[fraq_ind]-0.025
                                       }
                                       
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
                                     
                                     
                                     list(p_vals_D_fwer = p_vals_D_fwer, 
                                          p_vals_C_fwer = p_vals_C_fwer,
                                          p_vals_split_fwer = p_vals_split_fwer,
                                          p_vals_posi_fwer = p_vals_posi_fwer,
                                          fraq.Drysdale = fraq.vec.Drysdale[fraq_ind])
                                   }
  toc()
  stopCluster(cl)
  #Fetch p-values obtained from parallel computation
  results <- as.data.frame(results)
  p_vals_D_fwer <- results$p_vals_D_fwer
  p_vals_C_fwer <- results$p_vals_C_fwer
  p_vals_split_fwer <- results$p_vals_split_fwer
  p_vals_posi_fwer <- results$p_vals_posi_fwer
  #Compute confusion matrices, power and type1 error from parallel computation and average over all of them
  conf_matrices_D <- lapply(p_vals_D_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_C <- lapply(p_vals_C_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_split <-lapply(p_vals_split_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_posi <-lapply(p_vals_posi_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrix_all_D <- do.call(rbind, conf_matrices_D)
  conf_matrix_all_C <- do.call(rbind, conf_matrices_C)
  conf_matrix_all_split <- do.call(rbind, conf_matrices_split)
  conf_matrix_all_posi <- do.call(rbind, conf_matrices_posi)
  #Averaging over metrics
  test_res_D_avg <- colMeans(conf_matrix_all_D)
  test_res_C_avg <- colMeans(conf_matrix_all_C)
  test_res_split_avg <- colMeans(conf_matrix_all_split)
  test_res_posi_avg <- colMeans(conf_matrix_all_posi)
  power_avg_D <- test_res_D_avg[2]/(length(sel.index))
  power_avg_C <- test_res_C_avg[2]/(length(sel.index))
  power_avg_split <- test_res_split_avg[2]/(length(sel.index))
  power_avg_posi <- test_res_posi_avg[2]/(length(sel.index))
  type1_error_avg_D <- test_res_D_avg[1]/(p-length(sel.index))
  type1_error_avg_C <- test_res_C_avg[1]/(p-length(sel.index))
  type1_error_avg_split <- test_res_split_avg[1]/(p-length(sel.index))
  type1_error_avg_posi <- test_res_posi_avg[1]/(p-length(sel.index))
  #Calculate FWER
  H0T_Rej_any_D <- lapply(p_vals_D_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_C <- lapply(p_vals_C_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_split <- lapply(p_vals_split_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_posi <- lapply(p_vals_posi_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  FWER_D <- sum(do.call(rbind,H0T_Rej_any_D))/nsim
  FWER_C <- sum(do.call(rbind,H0T_Rej_any_C))/nsim
  FWER_split <- sum(do.call(rbind,H0T_Rej_any_split))/nsim
  FWER_posi <- sum(do.call(rbind,H0T_Rej_any_posi))/nsim
  #Printing results of one fraction to console
  cat("Results for fraction", fraq.vec[fraq_ind], ":\n")
  plot_conf_matrix(test_res_D_avg,"Drysdale's", nsim)
  plot_conf_matrix(test_res_C_avg,"Christoph's", nsim)
  cat("The average power of Drysdales p-values:", power_avg_D, "\n")
  cat("The average power of Christophs p-values:", power_avg_C,"\n")
  cat("The average power of splitting p-values:", power_avg_split,"\n")
  cat("The average power of posi p-values:", power_avg_posi,"\n")
  #cat("The average type 1 error of Drysdales p-values:", type1_error_avg_D,"\n")
  #cat("The average type 1 error of Christophs p-values:", type1_error_avg_C,"\n")
  cat("The FWER of Drysdales p-values:", FWER_D,"\n")
  cat("The FWER of Christophs p-values:", FWER_C,"\n")
  cat("The FWER of splitting p-values:", FWER_split,"\n")
  cat("The FWER of posi p-values:", FWER_posi,"\n")
  #Store everything to create plots later
  full_test_res_D[fraq_ind, ] <- test_res_D_avg
  full_test_res_C[fraq_ind, ] <- test_res_C_avg
  full_test_res_split[fraq_ind, ] <- test_res_split_avg
  full_test_res_posi[fraq_ind, ] <- test_res_posi_avg
  full_power_avg_D[fraq_ind] <- power_avg_D
  full_power_avg_C[fraq_ind] <- power_avg_C
  full_power_avg_split[fraq_ind] <- power_avg_split
  full_power_avg_posi[fraq_ind] <- power_avg_posi
  full_type1_error_avg_D[fraq_ind] <- type1_error_avg_D
  full_type1_error_avg_C[fraq_ind] <- type1_error_avg_C
  full_type1_error_avg_split[fraq_ind] <- type1_error_avg_split
  full_type1_error_avg_posi[fraq_ind] <- type1_error_avg_posi
  full_FWER_D[fraq_ind] <- FWER_D
  full_FWER_C[fraq_ind] <- FWER_C
  full_FWER_split[fraq_ind] <- FWER_split
  full_FWER_posi[fraq_ind] <- FWER_posi
  
  
  #Storing fraq.vec.Drysdale:
  avg_fraq.vec.Drysdale[fraq_ind] = do.call(sum, results$fraq.Drysdale)/nsim 
}

end.time <- Sys.time()
total.time <- end.time - start.time
cat("Total time needed for simulation:")
print(total.time)

#save.image(file='myEnvironment_nsim200_5active_sigma2_diff_fraqs_failed.RData')
#load('myEnvironment.RData')
# --------------- Create plots --------------

data_Power <- data.frame(
  Fraq=fraq.vec,
  "Avg Power Christoph" = full_power_avg_C,
  "Avg Power Drysdale" = full_power_avg_D,
  "Avg Power Splitting" = full_power_avg_split,
  "Avg Power Posi" = full_power_avg_posi
)

labels_fraqs_D<-c(rep("", 2), rep(avg_fraq.vec.Drysdale, 3))

data_Power_long <- tidyr::gather(data_Power, "Type", "Value", -Fraq)

PowerPlot<-ggplot(data_Power_long, aes(x = Fraq, y = Value, color = Type)) +
  geom_line() + geom_point() +
  geom_text(label=labels_fraqs_D, nudge_x=0.03, nudge_y=0.03) +
  labs(title = "Average Power",
       x = "Fraq", y = "Value") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5))
PowerPlot

ggsave("PowerPlot_diff_fraqs.png", plot = PowerPlot, width = 8, height = 6,
       units = "in", dpi = 300, bg = "#F0F0F0")


data_TypeI <- data.frame(
  Fraq=fraq.vec,
  "Avg Type I Error rate Christoph" = full_type1_error_avg_C,
  "Avg Type I Error rate Drysdale" = full_type1_error_avg_D,
  "Avg Type I Error rate Splitting" = full_type1_error_avg_split,
  "Avg Type I Error rate Posi" = full_type1_error_avg_posi
)

data_TypeI_long <- tidyr::gather(data_TypeI, "Type", "Value", -Fraq)

TypeIPlot<-ggplot(data_TypeI_long, aes(x = Fraq, y = Value, color = Type)) +
  geom_line() + geom_point() +
  geom_text(label=labels_fraqs_D, nudge_x=0.0002, nudge_y=0.0001) +
  labs(title = "Average Type I Error Rate",
       x = "Fraq", y = "Value") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5))

TypeIPlot

ggsave("TypeIPlot_diff_fraqs.png", plot = TypeIPlot, width = 8, height = 6,
       units = "in", dpi = 300, bg = "#F0F0F0")

#p_vals_posi_nofwer <- beta.posi(x, y, split=split, beta=beta_tmp,lambda=lambda, sigma=sigma_squ)$pvals_split

TypeIPlot
PowerPlot
