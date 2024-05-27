library(MASS)
library(glmnet)
library(tidyverse)
library(rsample)
library(ROCR)
library(truncnorm)

#setwd("~/R_sourcefiles/") #directory of sourcefiles
setwd("/Users/akira/Desktop/研究資料/RCVXMTL/git_code/R_sourcefiles/") #directory of sourcefiles

#load_sourcefiles
source("base_function.R")
source("estimator.R")
source("generate_data.R")
source("Initialize_R.R")


#settings
n_k <- 200 #number of samples in each task
n_kte <- 100 #number of samples for test data
n_ktr <- 50  #number of samples for training data
n_kv  <- n_k-(n_kte+n_ktr) #number of samples for validation data
p <- 100 #number of dimensions
T <- 150 #number of tasks
C <- 3  #number of true clusters
C_sd <- 10 #standard deviation for distribution of cluster center
w_var <- 1 #variance for cluster center 
ep_sd <- sqrt(5) #standard deviation of error term
rho_s <- 0 #correlation in design matrix

#settings concerning outlier parameters
o_var <- 5 #variance in outlier task
kappa <- 0.2  #probability of　occurence of outlier tasks
low_a <- 3 #lower bound for right truncated normal distribution 
up_b <- Inf #upper bound for right truncated normal distribution
trunc_mean <- 3 #mean of truncated normal distribution
rowwise_out = TRUE #TRUE means that all features are outlier for an outlier task. FALSE means that only non-zero features are outlier for an outlier task.

#data <- make_dataset_outlier_case1(T,n_k,p,C,C_sd,v_var,o_var,kappa,ep_sd,rho_s,low_a,up_b,trunc_mean ,rowwise_out=FALSE) #generating simulation data for Case1
data <- make_dataset_outlier_case2(T,n_k,p,C,C_sd,v_var,kappa,ep_sd,rho_s,rowwise_out=rowwise_out) #generating simulation data for Case2

X_list <-  data[[1]]
y_list <-  data[[2]]
W_ast_list <- data[[3]] #true regression coefficients
U_ast_list <- data[[4]] #true centroids
true_outlier_ind <- data[[7]] #index of true outlier tasks


W_ast <- matrix(unlist(W_ast_list), ncol = p, byrow = TRUE)

response = "gaussian"


#split samples into samples for learning and test
test_ind <- sample.int(n_k, size = n_kte, replace = FALSE , prob = NULL)

y_list_le <- y_list %>% purrr::map(~.x[-test_ind])
X_list_le <- X_list %>% purrr::map(~.x[-test_ind,])

y_list_te <- y_list %>% purrr::map(~.x[test_ind])
X_list_te <- X_list %>% purrr::map(~.x[test_ind,])

#split learning samples into samples for training and validation
train_ind <- sample.int(n_kv+n_ktr, size = n_ktr, replace = FALSE , prob = NULL)

y_list_tr <- y_list_le %>% purrr::map(~.x[train_ind])
X_list_tr <- X_list_le %>% purrr::map(~.x[train_ind,])

y_list_v <- y_list_le %>% purrr::map(~.x[-train_ind])
X_list_v <- X_list_le %>% purrr::map(~.x[-train_ind,])


#regularization term used for outlier parameters 
#thre_type <- g.soft.th #group_lasso
#thre_type <- gMCP #group_MCP
thre_type <- gSCAD #group_SCAD
#thre_type <- Hard_ridge #group_hard_ridge

#range_for_grid_search
lambda1_cand <- seq(0.2,0.4,by=0.2)
lambda2_cand <- seq(0.5,1,by=0.5)
lambda3_cand <- seq(2,3,by=1)

gamma = 3  #tuning parameter for non-convex penalties

Lambda <- expand.grid(lambda1_cand,lambda2_cand,lambda3_cand) #each row is a candidate of combinations about hyperparameters


W_ini <- t(mapply(Ini_glmnet, y_list_tr, X_list_tr,response=response,SIMPLIFY =T,alpha=1)) #initialization af regression coefficients by glmnet. alpha=1:lasso alpha=0:ridge, lambda is selected by CV.


#construct R 
K <- 5 # K-NN's K
R <- R_kNN(K,W_ini) #make R

#model selection
candidate_num <- nrow(Lambda) #number of model candidates
Verror_vec <- rep(1000,candidate_num) #vector that contains Validation error for each model candidate
for (j in 1:candidate_num)
{
  lambda1 <- Lambda[j,1]
  lambda2 <- Lambda[j,2]
  lambda3 <- Lambda[j,3]
  output <- MTLRRC_estimator(response=response,y_list_tr,X_list_tr,R,lambda1,lambda2,lambda3,gamma=gamma,thre_type = thre_type,rho=1,thre_ADMM=0.01,max_iter=200) #model estimation
  
  W_list <- output[[1]]
  if(response=="gaussian")
    Verror_vec[j] <- mean(mapply(cal_NMSE, y_list_v, X_list_v, W_list))
  if(response=="logistic")
    Verror_vec[j] <- mean(mapply(cal_Vlikelihood, y_list_v, X_list_v, W_list)) #calculating negative likelihood
}

#estimation by selected parameters
op_ind <- which.min(Verror_vec)
op_lambda1 <- Lambda[op_ind,1] #selected lambda1
op_lambda2 <- Lambda[op_ind,2] #selected lambda2
op_lambda3 <- Lambda[op_ind,3] #selected lambda3

output <- MTLRRC_estimator(response=response,y_list_tr,X_list_tr,R,op_lambda1,op_lambda2,op_lambda3,gamma=gamma,thre_type = thre_type,rho=1,thre_ADMM=0.01,max_iter=200) #model estimation with selected parameter

#estimated_model_parameters
W_hat_list <- output[[1]]
U_hat_list <- output[[2]]
O_hat_list <- output[[3]]


O_hat_mat <- do.call(rbind,O_hat_list)

est_outlier_ind <- NULL

for (m in 1:T)
{
  if( L2_norm(O_hat_list[[m]])>=0.01 )
  {
    est_outlier_ind <- append(est_outlier_ind,m)
  }
}


#calculate evaluation values
NMSE <- mean(mapply(cal_NMSE, y_list_te, X_list_te, W_hat_list))

W_RMSE <- sqrt(sum(mapply(sqL2_dif, W_hat_list,W_ast_list)))/T

U_hat_mat <- matrix(unlist(U_hat_list), ncol = p, byrow = TRUE)

TP <- length(intersect(est_outlier_ind,true_outlier_ind))
FP <- length(setdiff(est_outlier_ind,true_outlier_ind))
TPR <- TP/length(true_outlier_ind) # TP/(TP+FN)
FPR <- FP/length(setdiff(1:T,true_outlier_ind)) # FP/(FP+TN)




