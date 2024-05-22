# Clear all variables
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


source(hdi_adjustments_path)
source(carving_path)
source(sample_from_truncated_path)
source(tryCatchWE_path)
source(SNTN_distribution_path)
source(split_select_function_path)
source(carve_combined_path)



#-------------------- Toeplitz simulation from Multicarving paper ----------------------
n <- 100
p <- 200
rho <- 0.6
#fraq.vec <- c(0.5,0.6,0.7)
#fraq.vec <- c(0.7,0.8,0.9,0.95)
#fraq.vec <- c(0.5, 0.9)
fraq.vec <- c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99)
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
#sel.index <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)#active predictors
sel.index <- c(1,5,10,15,20)
beta <- rep(0, p)
beta[sel.index] <- 1
#RNGkind("Mersenne-Twister")#If we run multiple simulations in same R session, set this to true
set.seed(42) 
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta
#SNR = 4
SNR = 2
sigma_squ = drop(var(y.true)) / SNR

nsim <- 300
sig.level <- 0.05
new_fraq_threshold<-0
fraq.vec.comb<-fraq.vec

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
#Note posi fraq will be computed on same selected data as data splitting, whereas just
#posi data is computed on fraction = 1, so on all of the data
f <- length(fraq.vec)
full_test_res_D <- matrix(rep(0,4*f), nrow = f)
full_test_res_C <- matrix(rep(0,4*f), nrow = f)
full_test_res_split <- matrix(rep(0,4*f), nrow = f)
full_test_res_posi_fraq <- matrix(rep(0,4*f), nrow = f)
full_test_res_posi <- matrix(rep(0,4*f), nrow = f)
full_power_avg_D <- rep(0,f)
full_power_avg_C <- rep(0,f)
full_power_avg_split <- rep(0,f)
full_power_avg_posi_fraq <- rep(0,f)
full_power_avg_posi <- rep(0,f)
full_type1_error_avg_D <- rep(0,f)
full_type1_error_avg_C <- rep(0,f)
full_type1_error_avg_split <- rep(0,f)
full_type1_error_avg_posi_fraq <- rep(0,f)
full_type1_error_avg_posi <- rep(0,f)
full_FWER_D <- rep(0,f)
full_FWER_C <- rep(0,f)
full_FWER_split <- rep(0,f)
full_FWER_posi_fraq <- rep(0,f)
full_FWER_posi <- rep(0,f)
# To keep track of the average fraq used for beta^comb over all simulation runs:
avg_fraq.vec.comb<-rep(0,f)

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
     #Initialise fraq.vec.comb
     fraq.vec.comb <- fraq.vec
     
     
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
       p_vals_posi_fraq_fwer <- rep(1,p)
       print("0 variables where chosen by the lasso, but thats not a problem.t")
     }
     lambda <- split.select.list$lambda
     split <- split.select.list$split
     
     # --------- Extra variable selection for beta^comb ------------------
     
     split.select.list_D<-split.select.list
     beta_tmp_D<-split.select.list_D$beta
     lambda_D <- split.select.list_D$lambda
     split_D <- split.select.list_D$split
     
     # While the inverse can not be calculated: 
     while(sum(beta_tmp_D!=0)>n*(1-fraq.vec.comb[fraq_ind])){
       #Try new split with less observations for screening:
       
       if (counter_new_split>=new_fraq_threshold){
         fraq.vec.comb[fraq_ind]<-fraq.vec.comb[fraq_ind]-0.025
       }
       
       split.select.list_D <- split.select(x,y,fraction = fraq.vec.comb[fraq_ind])
       beta_tmp_D <- split.select.list_D$beta
       lambda_D <- split.select.list_D$lambda
       split_D <- split.select.list_D$split
       
       counter_new_split<-counter_new_split+1
     }
     
     if(sum(beta_tmp_D!=0)==0){
       empty_model_D <- TRUE
       p_vals_D_fwer <- rep(1,p)
       p_vals_split_fwer <- rep(1,p)
       p_vals_posi_fraq_fwer <- rep(1,p)
       print("0 variables where chosen by the lasso, but thats not a problem.")
     }
     #Compute pure p-values from combined carving estimator, carving estimator and regular splitting approach
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
     # I decided to include all of posi, split and combined carving in the case of whether
     # Or not beta^comb can be computed, since we work with fraq>0.7 anyways, so
     # beta^Posi not working shouldn't really be an issue
     
     if (!empty_model_D){
       carve_D <-carve.comb(x,y,split = split_D, beta = beta_tmp_D, lambda = lambda_D, sigma=sigma_squ)
       p_vals_D_nofwer <- carve_D$pvals
       p_vals_split_nofwer <- beta.split(x, y, split=split_D, beta = beta_tmp_D, sigma=sigma_squ)$pvals_split
       p_vals_posi_fraq_nofwer <- beta.posi(x, y, split=split_D, beta = beta_tmp_D,lambda=lambda_D, sigma=sigma_squ)$pvals#CHANGED LAMBDA TO LAMBDA_D
       
       chosen_D <- which(abs(beta_tmp_D)>0)
       model.size_D<- length(chosen_D)
       p_vals_D_fwer <- pmin(p_vals_D_nofwer * model.size_D, 1)
       p_vals_split_fwer <- pmin(p_vals_split_nofwer*model.size_D,1)
       p_vals_posi_fraq_fwer <- pmin(p_vals_posi_fraq_nofwer*model.size_D,1)
       
     }
     
     #new selection event on all of the data for posi100
     split.select.list.posi <- split.select(x,y,fraction = 1)
     beta_tmp_posi <- split.select.list.posi$beta
     lambda.posi <- split.select.list.posi$lambda
     split.posi <- split.select.list.posi$split
     chosen.posi <- which(abs(beta_tmp_posi)>0)
     model.size.posi <- length(chosen.posi)
     #handle empty model
     if(sum(beta_tmp_posi!=0)==0){
       p_vals_posi_fwer <- rep(1,p)
     }
     else{
       p_vals_posi_nofwer = beta.posi(x, y, split=split.posi, beta=beta_tmp_posi,lambda=lambda.posi, sigma=sigma_squ)$pvals
       p_vals_posi_fwer <- pmin(p_vals_posi_nofwer*model.size.posi,1)
       
     }
     
     
     list(p_vals_D_fwer = p_vals_D_fwer, 
          p_vals_C_fwer = p_vals_C_fwer,
          p_vals_split_fwer = p_vals_split_fwer,
          p_vals_posi_fraq_fwer = p_vals_posi_fraq_fwer,
          p_vals_posi_fwer = p_vals_posi_fwer,
          fraq.comb = fraq.vec.comb[fraq_ind])
   }
  toc()
  stopCluster(cl)
  #Fetch p-values obtained from parallel computation
  results <- as.data.frame(results)
  p_vals_D_fwer <- results$p_vals_D_fwer
  p_vals_C_fwer <- results$p_vals_C_fwer
  p_vals_split_fwer <- results$p_vals_split_fwer
  p_vals_posi_fraq_fwer <- results$p_vals_posi_fraq_fwer
  p_vals_posi_fwer <- results$p_vals_posi_fwer
  #Compute confusion matrices, power and type1 error from parallel computation and average over all of them
  conf_matrices_D <- lapply(p_vals_D_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_C <- lapply(p_vals_C_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_split <-lapply(p_vals_split_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_posi_fraq <-lapply(p_vals_posi_fraq_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_posi <-lapply(p_vals_posi_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrix_all_D <- do.call(rbind, conf_matrices_D)
  conf_matrix_all_C <- do.call(rbind, conf_matrices_C)
  conf_matrix_all_split <- do.call(rbind, conf_matrices_split)
  conf_matrix_all_posi_fraq <- do.call(rbind, conf_matrices_posi_fraq)
  conf_matrix_all_posi <- do.call(rbind, conf_matrices_posi)
  #Averaging over metrics
  test_res_D_avg <- colMeans(conf_matrix_all_D)
  test_res_C_avg <- colMeans(conf_matrix_all_C)
  test_res_split_avg <- colMeans(conf_matrix_all_split)
  test_res_posi_fraq_avg <- colMeans(conf_matrix_all_posi_fraq)
  test_res_posi_avg <- colMeans(conf_matrix_all_posi)
  power_avg_D <- test_res_D_avg[2]/(length(sel.index))
  power_avg_C <- test_res_C_avg[2]/(length(sel.index))
  power_avg_split <- test_res_split_avg[2]/(length(sel.index))
  power_avg_posi_fraq <- test_res_posi_fraq_avg[2]/(length(sel.index))
  power_avg_posi <- test_res_posi_avg[2]/(length(sel.index))
  type1_error_avg_D <- test_res_D_avg[1]/(p-length(sel.index))
  type1_error_avg_C <- test_res_C_avg[1]/(p-length(sel.index))
  type1_error_avg_split <- test_res_split_avg[1]/(p-length(sel.index))
  type1_error_avg_posi_fraq <- test_res_posi_fraq_avg[1]/(p-length(sel.index))
  type1_error_avg_posi <- test_res_posi_avg[1]/(p-length(sel.index))
  #Calculate FWER
  H0T_Rej_any_D <- lapply(p_vals_D_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_C <- lapply(p_vals_C_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_split <- lapply(p_vals_split_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_posi_fraq <- lapply(p_vals_posi_fraq_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_posi <- lapply(p_vals_posi_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  FWER_D <- sum(do.call(rbind,H0T_Rej_any_D))/nsim
  FWER_C <- sum(do.call(rbind,H0T_Rej_any_C))/nsim
  FWER_split <- sum(do.call(rbind,H0T_Rej_any_split))/nsim
  FWER_posi_fraq <- sum(do.call(rbind,H0T_Rej_any_posi_fraq))/nsim
  FWER_posi <- sum(do.call(rbind,H0T_Rej_any_posi))/nsim
  #Printing results of one fraction to console
  cat("Results for fraction", fraq.vec[fraq_ind], ":\n")
  plot_conf_matrix(test_res_D_avg,"Combined Carving", nsim)
  plot_conf_matrix(test_res_C_avg,"Carving", nsim)
  cat("The average power of combined carving p-values:", power_avg_D, "\n")
  cat("The average power of carving p-values:", power_avg_C,"\n")
  cat("The average power of splitting p-values:", power_avg_split,"\n")
  cat("The average power of posi p-values:", power_avg_posi,"\n")

  cat("The FWER of combined carvings p-values:", FWER_D,"\n")
  cat("The FWER of carvings p-values:", FWER_C,"\n")
  cat("The FWER of splitting p-values:", FWER_split,"\n")
  cat("The FWER of posi p-values:", FWER_posi,"\n")
  #Store everything to create plots later
  full_test_res_D[fraq_ind, ] <- test_res_D_avg
  full_test_res_C[fraq_ind, ] <- test_res_C_avg
  full_test_res_split[fraq_ind, ] <- test_res_split_avg
  full_test_res_posi_fraq[fraq_ind, ] <- test_res_posi_fraq_avg
  full_test_res_posi[fraq_ind, ] <- test_res_posi_avg
  full_power_avg_D[fraq_ind] <- power_avg_D
  full_power_avg_C[fraq_ind] <- power_avg_C
  full_power_avg_split[fraq_ind] <- power_avg_split
  full_power_avg_posi_fraq[fraq_ind] <- power_avg_posi_fraq
  full_power_avg_posi[fraq_ind] <- power_avg_posi
  full_type1_error_avg_D[fraq_ind] <- type1_error_avg_D
  full_type1_error_avg_C[fraq_ind] <- type1_error_avg_C
  full_type1_error_avg_split[fraq_ind] <- type1_error_avg_split
  full_type1_error_avg_posi_fraq[fraq_ind] <- type1_error_avg_posi_fraq
  full_type1_error_avg_posi[fraq_ind] <- type1_error_avg_posi
  full_FWER_D[fraq_ind] <- FWER_D
  full_FWER_C[fraq_ind] <- FWER_C
  full_FWER_split[fraq_ind] <- FWER_split
  full_FWER_posi_fraq[fraq_ind] <- FWER_posi_fraq
  full_FWER_posi[fraq_ind] <- FWER_posi
  
  
  #Storing fraq.vec.comb:
  avg_fraq.vec.comb[fraq_ind] = do.call(sum, results$fraq.comb)/nsim 
}
# Compute average results from PoSI estimator across all fractions, as it always worked with fraction 1
avg_power_posi <- mean(full_power_avg_posi)
avg_type1_error_posi <- mean(full_type1_error_avg_posi)
avg_fwer_posi <- mean(full_FWER_posi)

end.time <- Sys.time()
total.time <- end.time - start.time
cat("Total time needed for simulation:")
print(total.time)

save.image(file='Environment_diff_fraqs_s=5_SNR=2.RData')


#Need those NA's to integrate posi at fraction 1
data_Power <- data.frame(
  Fraq = c(fraq.vec, 1),
  "Carving" = c(full_power_avg_C, NA),
  "Combined Carving" = c(full_power_avg_D, NA),
  "Data Splitting" = c(full_power_avg_split, NA),
  "PoSI Fraq" = c(full_power_avg_posi_fraq, NA),
  "PoSI" = c(rep(NA, length(fraq.vec)), avg_power_posi)
)

FWER_points <- data.frame(
  Fraq = c(fraq.vec, 1),
  "Carving" = c(full_FWER_C, NA),
  "Combined Carving" = c(full_FWER_D, NA),
  "Data Splitting" = c(full_FWER_split, NA),
  "PoSI Fraq" = c(full_power_avg_posi_fraq, NA),
  "PoSI" = c(rep(NA, length(fraq.vec)), avg_fwer_posi)
)

# Convert data frames to long format
data_Power_long <- tidyr::gather(data_Power, "Type", "Value", -Fraq)
FWER_points_long <- tidyr::gather(FWER_points, "Type", "Value", -Fraq)

# Adjust fractions for points that overlap
FWER_points_adjusted <- FWER_points_long %>%
  group_by(Fraq, Value) %>%
  mutate(
    adjust_right = ifelse(duplicated(Value), 0.001, 0),
    adjust_left = ifelse(duplicated(Value, fromLast = TRUE), -0.001, 0),
    adjust_total = adjust_right + adjust_left,
    Fraq_adjusted = Fraq + adjust_total
  ) %>%
  ungroup()


PowerPlot <- ggplot(data_Power_long, aes(x = Fraq, y = Value, color = Type, linetype = Type, shape = Type), na.rm = TRUE) +
  geom_line(size = 1,na.rm = TRUE) +
  geom_hline(yintercept = sig.level, color = "black", linetype = "dashed") +
  geom_point(data = data_Power_long %>% filter(Fraq == 1), aes(x = 1, y = Value),size = 2.5, na.rm = TRUE) +
  geom_point(data = FWER_points_adjusted, aes(x = Fraq_adjusted, y = Value, color = Type, shape = Type), size = 2, na.rm = TRUE) +
  labs(title = "Average Power and FWER",
       x = "Fractions used for selection", y = "Value") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,0.8)) +
  scale_x_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "twodash", "blank"))+
  scale_shape_manual(values = c(0, 1, 2, 5, 7))+
  guides(linetype = guide_legend(override.aes = list(size = 3, alpha = 0.5)))


print(PowerPlot)
ggsave("diff_fraqs_s=5_SNR=2.png", plot = PowerPlot, width = 8, height = 6,
       units = "in", dpi = 300, bg = "#F0F0F0")

# Table for avg_fraq.vec.comb:
df_fraqs<-round(t(data.frame(original=fraq.vec, average=avg_fraq.vec.comb)),3)