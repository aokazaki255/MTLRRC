library(MASS)
library(glmnet)
library(rsample)
library(tidyverse)
library(ROCR)

#setwd("~/R_sourcefiles/") #directory of sourcefiles

#load_sourcefiles
source("base_function.R")
source("estimator.R")
source("generate_data.R")
source("Initialize_R.R")


dataset = "landmine"
#dataset = "school"

if(dataset == "school")
{
  T <- 139 #number of tasks
  response <- "gaussian"
  dataset_path <- "" #set your dataset directory
}

if(dataset=="landmine")
{
  T <-29 
  response <-"binomial"
  dataset_path <- ""  #set your dataset directory
}


tr_prop <- 0.6 #proportion of train_data
v_prop  <- 0.2 #proportion of validation_data
te_prop <- 1-(tr_prop+v_prop) #proportion of test_data
data_list <- make_real_dataset(tr_prop,v_prop,te_prop,data=dataset ,path = dataset_path)

tr_data <- data_list[[1]] 
v_data <- data_list[[2]]
te_data <- data_list[[3]]

y_list_v <-  v_data[[1]]
X_list_v <- v_data[[2]]
y_list_tr <-  tr_data[[1]]
X_list_tr <- tr_data[[2]]


#regularization term used for outlier parameters 
#thre_type <- g.soft.th #group_lasso
#thre_type <- gMCP #group_MCP
thre_type <- gSCAD #group_SCAD
#thre_type <- Hard_ridge #group_hard_ridge

#range_for_grid_search
lambda1_cand <- seq(0.2,0.4,by=0.2)
lambda2_cand <- seq(0.5,1,by=0.5)
lambda3_cand <- seq(0.5,1,by=0.5)
gamma_cand <- seq(2.01,2.51, by = 0.5)


Lambda <- expand.grid(lambda1_cand,lambda2_cand,lambda3_cand,gamma_cand) #each row is a candidate of combinations about hyperparameters


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
  gamma <- Lambda[j,4]
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
op_gamma <- Lambda[op_ind,4] #selected gamma

output <- MTLRRC_estimator(response=response,y_list_tr,X_list_tr,R,op_lambda1,op_lambda2,op_lambda3,gamma=op_gamma,thre_type = thre_type,rho=1,thre_ADMM=0.01,max_iter=200) #model estimation with selected parameter

#estimated_model_parameters
W_hat_list <- output[[1]]
U_hat_list <- output[[2]]
O_hat_list <- output[[3]]



y_list_te <- te_data[[1]]
X_list_te <- te_data[[2]]

y_list_hat <- vector("list",length = T)
n_list_te <- lapply(y_list_te, length)


#calculate evaluation value
p <- ncol(X_list_te[[1]])
if(response=="gaussian")
{
  NMSE <- mean(mapply(cal_NMSE, y_list_te, X_list_te, W_hat_list)) #normalized mean squared error
  W_hat_mat <- matrix(unlist(W_hat_list),ncol=p,byrow=TRUE)
  print(NMSE)
}
  
if(response=="binomial")
{
  AUC_vec <- rep(0,T)
  for(k in 1:T)
  {
    n <- length(y_list_te[[k]])
    y_hat_te <- rep(0,n)
    for (i in 1:n)
    {
      y_hat_te[i] <- pi(c(1,X_list_te[[k]][i,]),W_hat_list[[k]])
    }
    pred <- prediction(y_hat_te, y_list_te[[k]])
    auc.tmp <- performance(pred,"auc")
    auc <- as.numeric(auc.tmp@y.values)
    AUC_vec[k] <- auc
  }
  AUC <- mean(AUC_vec)
  W_hat_mat <- matrix(unlist(W_hat_list),ncol=p+1,byrow=TRUE)
  print(AUC)
}

U_hat_mat <- matrix(unlist(U_hat_list),ncol=p,byrow=TRUE)
O_hat_mat <- matrix(unlist(O_hat_list),ncol=p,byrow=TRUE)



