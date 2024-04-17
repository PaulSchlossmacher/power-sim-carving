# ----------- Dependencies

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
carve_linear_clean_path<-paste(Local_path, "/carve_linear_clean.R", sep="")

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
source(carve_linear_clean_path)

#--------------------------


#-------------------- Toeplitz Carving simulation from Christoph ----------------------
n <- 100
p <- 200
rho <- 0.6
fraq = 0.7
#toeplitz takes the first column of the desired toeplitz design and creates the whole function, here a sequence from 0 to p-1
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
sel.index <- c(1, 5, 10, 15, 20)#active predictors
beta_0 <- rep(0, p)#initialize beta as all zeros
beta_0[sel.index] <- 1#put ones at active predictor positions
sparsity <- 5
set.seed(41) # to make different methods comparable, fix the x-matrix
x <- mvrnorm(n, rep(0, p), Cov)#sample X from multivariate normal distribution
y.true <- x %*% beta_0
SNR <- 1.713766 # value created for Toeplitz 0.6
sigma <- 2
y <- y.true + sigma * rnorm(n)


split.select.list <- split.select(x,y,fraction = fraq)
lambda <- split.select.list$lambda
split <- split.select.list$split
beta <- split.select.list$beta


#Normalize x and y before starting:
y<-(y-mean(y))
for (j in 1:dim(x)[2]){
  xjbar<-mean(x[,j])
  sigma_j_squ<-sum((x[,j]-xjbar)^2)/(length(x[,j])-1)
  for (i in 1:dim(x)[1]){
    x[i,j]<-(x[i,j]-xjbar)/sqrt(sigma_j_squ)
  }
}

#Split the data
n <- length(y)
p <- length(beta)
n.a <- length(split)
n.b <- n-n.a
x.a <- x[split, ]
y.a <- y[split, ]
x.b <- x[-split, ]
y.b <- y[-split, ]
sigma_squ <- sigma
#Sigma gets chosen in accordance with the distribution of y_A in Lemma 3.2
Sigma <- diag(n.a)*sigma_squ
c1<-n.b/n
c2<-n.a/n

chosen <-  which(abs(beta) > 0) # selected variables
s <- length(chosen)
b.signs <- sign(beta[chosen])


if (s == 0)
  stop("0 variables were chosen by the Lasso")

#extract active variables from both splits
x.Ma <- x.a[, chosen]#(n.a x s)
x.Mb <- x.b[, chosen]#(n.b x s)

#extract inactive variables from both splits
x_Ma <- x.a[, -chosen]#(n.a x (p-s))
x_Mb <- x.b[, -chosen]

#compute the Moore-Penrose inverse of active variables in both splits
x.Ma.i <- ginv(x.Ma)#(s x na)
x.Ma.ti <- ginv(t(x.Ma))
x.Mb.i <- ginv(x.Mb)

#compute the projection matrix onto active columns
p.Ma <- x.Ma %*% x.Ma.i#(n.a x n.a)

#compute the beta_carve from drysdales paper
beta_carve_D <- ((n.a/n)*x.Ma.i %*% y.a + (n.b/n)*x.Mb.i %*% y.b)

#Get inactive affine constraints on split A, as this is the group where we are doing POSI(Lee et al. page 8)
A.0.up <- (t(x_Ma) %*% (diag(n.a) - p.Ma))#(p-s) x n.a
A.0 <- 1/lambda*rbind(A.0.up, -A.0.up) #2*(p-s) x n.a
b.0.temp <-t(x_Ma) %*% x.Ma.ti %*% b.signs
b.0.up <- rep(1,p-s) - b.0.temp
b.0.lo <- rep(1,p-s) + b.0.temp
b.0 <- rbind(b.0.up, b.0.lo)#2*(p-s) x n.a

#Get active affine constraints on split A
C <- solve(crossprod(x.Ma,x.Ma))#(X_Ma^T*X_Ma)^-1
A.1 <- -diag(x = b.signs, nrow = s) %*% x.Ma.i # (s x n.a)
b.1 <- -lambda*diag(x = b.signs, nrow = s) %*% C %*% b.signs#Here Christoph differentiates for intercept

A <- rbind(A.0, A.1)# (2p-s) x n.a
b <- rbind(b.0,b.1)

#Following a mix of Lee (https://github.com/selective-inference/R-software/blob/master/selectiveInference/R/funs.fixed.R) from line 231
#and Drysdale (https://github.com/ErikinBC/sntn/blob/main/sntn/_lasso.py) from line 195
vup <- rep(0,s)
vlo <- rep(0,s)
v0 <- rep(0,s)

keep_zero_ind <- c()
keep_zero_resid <- list()
keep_zero_den <- list()
keep_zero_A <- list()
keep_zero_c <- list()

keep_ind_vlo <- list()
norm_consts <- rep(0,s)
for (i in 1:s){
  v.i <- x.Ma.i[i,]
  v.i.norm <- sqrt(sum(v.i^2))
  
  if(TRUE){
    eta <- b.signs[i]*v.i/v.i.norm
  }
  else {
    eta <- x.Ma.i[i,]
  }
  c <- (Sigma %*% eta) / as.numeric((t(eta) %*% Sigma %*% eta))
  z <- (diag(n.a) - c %*% t(eta)) %*% y.a
  den <- A%*%c
  resid <- b-A %*% z
  #We do not consider the set V^0(z) as defined in Lee p.10, because
  #Drysdale does not do so either
  ind.vup <- (den > 0)
  ind.vup <- which(ind.vup == TRUE)[which(resid[ind.vup]>0)]
  ind.vlo <- (den < 0)
  ind.vlo <- which(ind.vlo == TRUE)[which(resid[ind.vlo]>0)]
  ind.v0 <- (den == 0)
 
  if (any(ind.vup)){
    vup[i] <- min(resid[ind.vup]/den[ind.vup])
  }else {
    vup[i] <- Inf
  }
  
  if (any(ind.vlo)){
    vlo[i] <- max(resid[ind.vlo]/den[ind.vlo])
    keep_ind_vlo[[i]] = ind.vlo
  }
  else {
    vlo[i] <- -Inf
  }
  
  
  # ------- For testing V0(z)
  #According to Lemma 5.1 in Lee's paper, {Ay<=b}={V-(z)<=eta^T y <=V+(z), V0(z)>=0}
  
  if (any(ind.v0)){
    v0[i] <- min(resid[ind.v0])
    keep_zero_ind <- c(keep_zero_ind, i)
    resid_zero_name <- paste0("resid",i)
    den_zero_name <- paste0("den",i)
    A_zero_name <- paste0("A",i)
    c_zero_name <- paste0("c",i)
    keep_zero_resid[[resid_zero_name]] <- resid
    keep_zero_den[[den_zero_name]] <- den
    keep_zero_A[[A_zero_name]] <- A 
    keep_zero_c[[c_zero_name]] <- c
  }else {
    
    # I guess the choice here could be disputed, but 0 makes the most sense to me
    v0[i] <- 0
  }
  
  norm_consts[i] <- v.i.norm
}

v0
vlo
vup
# for (i in 1:length(keep_zero_den)){
#   print(which(keep_zero_den[[i]] == 0))
# }
# 
# results <- sapply(keep_zero_den, function(x) {
#   # Get indices where elements are 0
#   zero_indices <- which(x == 0)
#   
#   # Check if any of these indices are greater than the number of rows in A.0
#   any(zero_indices > dim(A.0)[1])
# })

