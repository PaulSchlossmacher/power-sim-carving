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
#fraq.vec <- c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99)
fraq.vec <- c(0.7,0.75,0.8,0.85,0.9,0.95,0.99)
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
#sel.index <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)#active predictors
sel.index <- c(1,5,10,15,20)
sparsity <- length(sel.index)
beta <- rep(0, p)
beta[sel.index] <- 1
#RNGkind("Mersenne-Twister")#If we run multiple simulations in same R session, set this to true
set.seed(42) 
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta
#SNR = 4
SNR <- 2
sigma_squ <- drop(var(y.true)) / SNR
sigma <- sqrt(sigma_squ)
#sigma_squ = 1
#sigma_squ <- 2 #variance used in some other simulations without fixing SNR
nsim <- 50
sig.level <- 0.05
flexible_selection  <- FALSE #should we allow more selections for data splitting approach
flexible_selection_count <- 5


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
full_test_res_posi <- matrix(rep(0,4*f), nrow = f)
full_test_res_posi_fraq <- matrix(rep(0,4*f), nrow = f)
full_test_res_sat <- matrix(rep(0,4*f), nrow = f)
full_power_avg_D <- rep(0,f)
full_power_avg_C <- rep(0,f)
full_power_avg_split <- rep(0,f)
full_power_avg_posi <- rep(0,f)
full_power_avg_posi_fraq <- rep(0,f)
full_power_avg_sat <- rep(0,f)
full_type1_error_avg_D <- rep(0,f)
full_type1_error_avg_C <- rep(0,f)
full_type1_error_avg_split <- rep(0,f)
full_type1_error_avg_posi <- rep(0,f)
full_type1_error_avg_posi_fraq <- rep(0,f)
full_type1_error_avg_sat <- rep(0,f)
full_FWER_D <- rep(0,f)
full_FWER_C <- rep(0,f)
full_FWER_split <- rep(0,f)
full_FWER_posi <- rep(0,f)
full_FWER_posi_fraq <- rep(0,f)
full_FWER_sat <- rep(0,f)

log.vec <- vector("list",length(fraq.vec))
#Should count number of repeated selection events to make combined carving estimator work over nsim rounds and a given fraction
repeated_selections <- rep(0,f)
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
    #get different selection events
    logs <- c()
    select.again <- TRUE
    empty_model <- FALSE
    splitting_estimator_failed <- FALSE
    flexible_selection_par <- flexible_selection#Should be reset before every round
    select.again.counter = 0
    while(select.again){
      if (select.again.counter >= flexible_selection_count){
        warning("Tried to many selection events and not one of them was conformable for beta_comb and beta_split")
        flexible_selection_par <- FALSE
      }
      select.again <- FALSE
      y <- y.true + sqrt(sigma_squ) * rnorm(n)
      #Normalize y:
      y<-(y-mean(y))
      split.select.list <- split.select(x,y,fraction = fraq.vec[fraq_ind])
      beta_tmp <- split.select.list$beta
      if(sum(beta_tmp!=0)==0){
        empty_model <- TRUE
        p_vals_D_fwer <- rep(1,p)
        p_vals_C_fwer <- rep(1,p)
        p_vals_posi_fraq_fwer <- rep(1,p)
        p_vals_split_fwer <- rep(1,p)
        p_vals_sat_fwer <- rep(1,p)
        logs <- c(logs,"0 variables where chosen by the lasso, but thats not a problem.")
      }
      lambda <- split.select.list$lambda
      split <- split.select.list$split
      if(sum(beta_tmp!=0)>min(n*fraq.vec[fraq_ind], n*(1-fraq.vec[fraq_ind])) && flexible_selection_par){
        select.again <- TRUE
        select.again.counter <- select.again.counter + 1
        logs <- c(logs,"Need to split again because we selected more variables than beta_D & beta_split can handle")
      }
      else if(sum(beta_tmp!=0)>min(n*fraq.vec[fraq_ind], n*(1-fraq.vec[fraq_ind])) && !flexible_selection_par){
        p_vals_D_fwer <- rep(1,p)
        p_vals_split_fwer <- rep(1,p)
        p_vals_posi_fraq_fwer <- rep(1,p)
        splitting_estimator_failed <- TRUE
        logs <- c(logs,"Set splitting and carve.comb p_vals to 1, because we selected more variables than beta_D & beta_split can handle")
        
      }
    }
    
    
    if(!empty_model && !splitting_estimator_failed){
      #Compute pure p-values from combined carving estimator, carving estimator and regular splitting approach
      carve_D <-carve.comb(x,y,split = split, beta = beta_tmp, lambda = lambda, sigma_squ=sigma_squ)
      p_vals_D_nofwer <- carve_D$pvals
      
      carve_C <- carve.lasso(X = x, y = y, ind = split, beta = beta_tmp, tol.beta = 0, sigma = sigma,
                             lambda = lambda,FWER = FALSE, intercept = FALSE,selected=TRUE, verbose = FALSE)
      p_vals_C_nofwer<-carve_C$pv
      
      #Note: selected=FALSE for saturated model
      carve_sat <-carve.lasso(X = x, y = y, ind = split, beta = beta_tmp, tol.beta = 0, sigma = sigma,
                              lambda = lambda,FWER = FALSE, intercept = FALSE, selected=FALSE, verbose = FALSE)
      p_vals_sat_nofwer<-carve_sat$pv
      
      p_vals_split_nofwer <- beta.split(x, y, split=split, beta=beta_tmp, sigma_squ=sigma_squ)$pvals_split
      
      p_vals_posi_fraq_nofwer <- beta.posi(x, y, split=split, beta=beta_tmp,lambda=lambda, sigma_squ=sigma_squ)$pvals
      
      #carve_C only returns the p-values of the coefficients determined by the selection event, hence we assign them at the appropriate positions
      p_vals_comp_C<-rep(1,p)
      chosen <- which(abs(beta_tmp)>0)
      p_vals_comp_C[chosen] <- p_vals_C_nofwer
      p_vals_comp_sat<-rep(1,p)
      p_vals_comp_sat[chosen] <- p_vals_sat_nofwer
      
      #Add FWER control with Bonferroni correction
      model.size <- length(chosen)
      p_vals_D_fwer <- pmin(p_vals_D_nofwer * model.size, 1)
      p_vals_C_fwer <- pmin(p_vals_comp_C * model.size, 1)
      p_vals_split_fwer <- pmin(p_vals_split_nofwer*model.size,1)
      p_vals_posi_fraq_fwer <- pmin(p_vals_posi_fraq_nofwer*model.size,1)
      p_vals_sat_fwer <- pmin(p_vals_comp_sat*model.size,1)
      
    }
    #we selected something but data splitting failed, so we compute only the estimators that dont depend on it
    else if(!empty_model){
      carve_C <- carve.lasso(X = x, y = y, ind = split, beta = beta_tmp, tol.beta = 0, sigma = sigma,
                             lambda = lambda,FWER = FALSE, intercept = FALSE,selected=TRUE, verbose = FALSE)
      p_vals_C_nofwer<-carve_C$pv
      
      carve_sat <-carve.lasso(X = x, y = y, ind = split, beta = beta_tmp, tol.beta = 0, sigma = sigma,
                              lambda = lambda,FWER = FALSE, intercept = FALSE, selected=FALSE, verbose = FALSE)
      p_vals_sat_nofwer<-carve_sat$pv

      p_vals_comp_C<-rep(1,p)
      chosen <- which(abs(beta_tmp)>0)
      p_vals_comp_C[chosen] <- p_vals_C_nofwer
      
      p_vals_comp_sat<-rep(1,p)
      p_vals_comp_sat[chosen] <- p_vals_sat_nofwer
      
      #Add FWER control with Bonferroni correction
      model.size <- length(chosen)
      p_vals_C_fwer <- pmin(p_vals_comp_C * model.size, 1)
      p_vals_sat_fwer <- pmin(p_vals_comp_sat*model.size,1)
    }
    
    #new selection event on all of the data for posi
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
      p_vals_posi_nofwer = beta.posi(x, y, split=split.posi, beta=beta_tmp_posi,lambda=lambda.posi, sigma_squ=sigma_squ)$pvals
      p_vals_posi_fwer <- pmin(p_vals_posi_nofwer*model.size.posi,1)
      
    }
    
    list(p_vals_D_fwer = p_vals_D_fwer, 
         p_vals_C_fwer = p_vals_C_fwer,
         p_vals_sat_fwer = p_vals_sat_fwer,
         p_vals_split_fwer = p_vals_split_fwer,
         p_vals_posi_fwer = p_vals_posi_fwer,
         p_vals_posi_fraq_fwer = p_vals_posi_fraq_fwer,
         select.again.counter = select.again.counter,
         logs = logs)
    
  }
  toc()
  stopCluster(cl)
  
  #Fetch p-values obtained from parallel computation
  results <- as.data.frame(results)
  p_vals_D_fwer <- results$p_vals_D_fwer
  p_vals_C_fwer <- results$p_vals_C_fwer
  p_vals_split_fwer <- results$p_vals_split_fwer
  p_vals_posi_fwer <- results$p_vals_posi_fwer
  p_vals_posi_fraq_fwer <- results$p_vals_posi_fraq_fwer
  p_vals_sat_fwer <- results$p_vals_sat_fwer
  select.again.counter <- do.call(rbind,results$select.again.counter)
  
  #Compute confusion matrices, power and type1 error from parallel computation and average over all of them
  conf_matrices_D <- lapply(p_vals_D_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_C <- lapply(p_vals_C_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_split <-lapply(p_vals_split_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_posi <-lapply(p_vals_posi_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_posi_fraq <-lapply(p_vals_posi_fraq_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrices_sat <-lapply(p_vals_sat_fwer, function(p_vals) conf_matrix(p_vals, sig.level = sig.level, beta = beta))
  conf_matrix_all_D <- do.call(rbind, conf_matrices_D)
  conf_matrix_all_C <- do.call(rbind, conf_matrices_C)
  conf_matrix_all_split <- do.call(rbind, conf_matrices_split)
  conf_matrix_all_posi <- do.call(rbind, conf_matrices_posi)
  conf_matrix_all_posi_fraq <- do.call(rbind, conf_matrices_posi_fraq)
  conf_matrix_all_sat <- do.call(rbind, conf_matrices_sat)
  
  #Averaging over metrics
  test_res_D_avg <- colMeans(conf_matrix_all_D)
  test_res_C_avg <- colMeans(conf_matrix_all_C)
  test_res_split_avg <- colMeans(conf_matrix_all_split)
  test_res_posi_avg <- colMeans(conf_matrix_all_posi)
  test_res_posi_fraq_avg <- colMeans(conf_matrix_all_posi_fraq)
  test_res_sat_avg <- colMeans(conf_matrix_all_sat)
  power_avg_D <- test_res_D_avg[2]/(length(sel.index))
  power_avg_C <- test_res_C_avg[2]/(length(sel.index))
  power_avg_split <- test_res_split_avg[2]/(length(sel.index))
  power_avg_posi <- test_res_posi_avg[2]/(length(sel.index))
  power_avg_posi_fraq <- test_res_posi_fraq_avg[2]/(length(sel.index))
  power_avg_sat <- test_res_sat_avg[2]/(length(sel.index))
  type1_error_avg_D <- test_res_D_avg[1]/(p-length(sel.index))
  type1_error_avg_C <- test_res_C_avg[1]/(p-length(sel.index))
  type1_error_avg_split <- test_res_split_avg[1]/(p-length(sel.index))
  type1_error_avg_posi <- test_res_posi_avg[1]/(p-length(sel.index))
  type1_error_avg_posi_fraq <- test_res_posi_fraq_avg[1]/(p-length(sel.index))
  type1_error_avg_sat <- test_res_sat_avg[1]/(p-length(sel.index))
  
  #Calculate FWER
  H0T_Rej_any_D <- lapply(p_vals_D_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_C <- lapply(p_vals_C_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_split <- lapply(p_vals_split_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_posi <- lapply(p_vals_posi_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_posi_fraq <- lapply(p_vals_posi_fraq_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  H0T_Rej_any_sat <- lapply(p_vals_sat_fwer, function(p_vals) any(p_vals<=sig.level & beta==0))
  FWER_D <- sum(do.call(rbind,H0T_Rej_any_D))/nsim
  FWER_C <- sum(do.call(rbind,H0T_Rej_any_C))/nsim
  FWER_split <- sum(do.call(rbind,H0T_Rej_any_split))/nsim
  FWER_posi <- sum(do.call(rbind,H0T_Rej_any_posi))/nsim
  FWER_posi_fraq <- sum(do.call(rbind,H0T_Rej_any_posi_fraq))/nsim
  FWER_sat <- sum(do.call(rbind,H0T_Rej_any_sat))/nsim
  
  #Printing results of one fraction to console for carving,combined carving and data splitting
  cat("Results for fraction", fraq.vec[fraq_ind], ":\n")
  plot_conf_matrix(test_res_D_avg,"Combined Carving", nsim)
  plot_conf_matrix(test_res_C_avg,"Carving", nsim)
  cat("The average power of combined carving p-values:", power_avg_D, "\n")
  cat("The average power of carving p-values:", power_avg_C,"\n")
  cat("The average power of splitting p-values:", power_avg_split,"\n")
  cat("The FWER of combined carvings p-values:", FWER_D,"\n")
  cat("The FWER of carvings p-values:", FWER_C,"\n")
  cat("The FWER of splitting p-values:", FWER_split,"\n")

  #Store everything to create plots later
  full_test_res_D[fraq_ind, ] <- test_res_D_avg
  full_test_res_C[fraq_ind, ] <- test_res_C_avg
  full_test_res_split[fraq_ind, ] <- test_res_split_avg
  full_test_res_posi[fraq_ind, ] <- test_res_posi_avg
  full_test_res_posi_fraq[fraq_ind, ] <- test_res_posi_fraq_avg
  full_test_res_sat[fraq_ind, ] <- test_res_sat_avg
  full_power_avg_D[fraq_ind] <- power_avg_D
  full_power_avg_C[fraq_ind] <- power_avg_C
  full_power_avg_split[fraq_ind] <- power_avg_split
  full_power_avg_posi[fraq_ind] <- power_avg_posi
  full_power_avg_posi_fraq[fraq_ind] <- power_avg_posi_fraq
  full_power_avg_sat[fraq_ind] <- power_avg_sat
  full_type1_error_avg_D[fraq_ind] <- type1_error_avg_D
  full_type1_error_avg_C[fraq_ind] <- type1_error_avg_C
  full_type1_error_avg_split[fraq_ind] <- type1_error_avg_split
  full_type1_error_avg_posi[fraq_ind] <- type1_error_avg_posi
  full_type1_error_avg_posi_fraq[fraq_ind] <- type1_error_avg_posi_fraq
  full_type1_error_avg_sat[fraq_ind] <- type1_error_avg_sat
  full_FWER_D[fraq_ind] <- FWER_D
  full_FWER_C[fraq_ind] <- FWER_C
  full_FWER_split[fraq_ind] <- FWER_split
  full_FWER_posi[fraq_ind] <- FWER_posi
  full_FWER_posi_fraq[fraq_ind] <- FWER_posi_fraq
  full_FWER_sat[fraq_ind] <- FWER_sat
  repeated_selections[fraq_ind] <- sum(select.again.counter)
  log.vec[[fraq_ind]] <- unlist(results$logs)
}
# Compute average results from PoSI estimator across all fractions, as it always worked with fraction 1
avg_power_posi <- mean(full_power_avg_posi)
avg_type1_error_posi <- mean(full_type1_error_avg_posi)
avg_fwer_posi <- mean(full_FWER_posi)

end.time <- Sys.time()
total.time <- end.time - start.time
cat("Simulation was successful")
print(total.time)
#Print information about repeated selections
rep_select_df <- data.frame(fractions = fraq.vec, repeated_selections = repeated_selections)
print(rep_select_df)

#Save or load existing simulation environments
if (flexible_selection){
  filename <- sprintf("m=%d_SNR=%.1f_allowed_fails", sparsity, SNR)
}else{
  filename <- sprintf("m=%d_SNR=%.1f", sparsity, SNR)
}
filename <- gsub(".0", "", filename)
filename <- sub("\\.", ",", filename)
environment_name <- paste0("Environment_",filename,".RData")
plot_name <- paste0("main_plot_",filename,".png")
save.image(file=environment_name)
#load("simulation_environments/Environment_s=5_SNR=2.RData")


# --------------- Create main power plots --------------
if (flexible_selection){
  main_title <- sprintf("m = %d, SNR = %.1f, max_rep_sel = %d", sparsity, SNR, flexible_selection_count)
}else{
  main_title <- sprintf("m = %d, SNR = %.1f", sparsity, SNR)
}
main_title <- gsub(".0", "", main_title)
#Need those NA's to integrate posi at fraction 1
data_Power <- data.frame(
  Fraq = c(fraq.vec, 1),
  "Carving" = c(full_power_avg_C, NA),
  "Carving Sat." = c(full_power_avg_sat,NA),
  "Combined Carving" = c(full_power_avg_D, NA),
  "Data Splitting" = c(full_power_avg_split, NA),
  "PoSI" = c(rep(NA, length(fraq.vec)), avg_power_posi)
)

FWER_points <- data.frame(
  Fraq = c(fraq.vec, 1),
  "Carving" = c(full_FWER_C, NA),
  "Carving Sat." = c(full_FWER_sat, NA),
  "Combined Carving" = c(full_FWER_D, NA),
  "Data Splitting" = c(full_FWER_split, NA),
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
  geom_point(data = data_Power_long %>% filter(Fraq == 1), aes(x = 1, y = Value),size = 2.5,stroke = 1.2, na.rm = TRUE) +
  geom_point(data = FWER_points_adjusted, aes(x = Fraq_adjusted, y = Value, color = Type, shape = Type), size = 2.5, stroke = 1.2,alpha = 0.7, na.rm = TRUE) +
  labs(title = main_title,
       x = "Fractions used for selection", y = "Average power(-) and FWER(o)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    #legend.position = "none",
    legend.position = c(0.85,0.85),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.key.size = unit(2, "lines"),
    legend.background = element_rect( color = "black"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm")
  )+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "twodash", "blank"))+
  scale_shape_manual(values = c(0, 1, 2, 5, 7))+
  guides(linetype = guide_legend(override.aes = list(size = 3,alpha = 0.5)))


print(PowerPlot)
ggsave(plot_name, plot = PowerPlot, width = 8, height = 6,
       units = "in", dpi = 300, bg = "#F0F0F0")


# --------------- Create plots for visualization of power composition in combined carving estimator --------------

# -------- Create plots for visualization of power composition in combined carving estimator ---------
#plot_name2 <- paste0("combined_plot_",filename,".png")

#Need those NA's to integrate posi at fraction 1
# data_Power2 <- data.frame(
#   Fraq = fraq.vec,
#   "Combined Carving" = full_power_avg_D,
#   "Data Splitting" = full_power_avg_split,
#   "PoSI" = full_power_avg_posi_fraq
#   )
# 
# FWER_points2 <- data.frame(
#   Fraq = fraq.vec,
#   "Combined Carving" = full_FWER_D,
#   "Data Splitting" = full_FWER_split,
#   "PoSI" = full_FWER_posi_fraq
# )
# 
# # Convert data frames to long format
# data_Power2_long <- tidyr::gather(data_Power2, "Type", "Value", -Fraq)
# FWER_points2_long <- tidyr::gather(FWER_points2, "Type", "Value", -Fraq)
# 
# # Adjust fractions for points that overlap
# FWER_points_adjusted2 <- FWER_points2_long %>%
#   group_by(Fraq, Value) %>%
#   mutate(
#     adjust_right = ifelse(duplicated(Value), 0.001, 0),
#     adjust_left = ifelse(duplicated(Value, fromLast = TRUE), -0.001, 0),
#     adjust_total = adjust_right + adjust_left,
#     Fraq_adjusted = Fraq + adjust_total
#   ) %>%
#   ungroup()
# 
# PowerPlot2 <- ggplot(data_Power2_long, aes(x = Fraq, y = Value, color = Type, linetype = Type, shape = Type)) +
#   geom_line(size = 1) +
#   geom_hline(yintercept = sig.level, color = "black", linetype = "dashed") +
#   geom_point(data = FWER_points_adjusted2, aes(x = Fraq_adjusted, y = Value, color = Type, shape = Type),  size = 2.5, stroke = 1.2,alpha = 0.7) +
#   labs(title = "m = 5, SNR = 2",
#        x = "Fractions used for selection", y = "Average power(-) and FWER(o)") +
#   theme_minimal() + 
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 25),
#     axis.title.x = element_text(size = 20),
#     axis.title.y = element_text(size = 20),
#     #legend.position = "none",
#     legend.position = c(0.85,0.85),
#     legend.text = element_text(size = 14),
#     legend.title = element_blank(),
#     legend.key.size = unit(2, "lines"),
#     legend.background = element_rect( color = "black"),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     axis.line = element_line(color = "black"),
#     axis.ticks = element_line(color = "black"),
#     axis.ticks.length = unit(0.25, "cm")
#   )+
#   scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
#   scale_x_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
#   scale_linetype_manual(values = c("solid","dotdash", "longdash"))+
#   scale_shape_manual(values = c(2, 3, 1))+
#   guides(linetype = guide_legend(override.aes = list(size = 3,alpha = 0.5)))
# 
# print(PowerPlot2)
# ggsave(plot_name2, plot = PowerPlot2, width = 8, height = 6,
#        units = "in", dpi = 300, bg = "#F0F0F0")



