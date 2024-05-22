#calculate MSE
sqL2_dif <- function(x,y)
{
  return(sqL2_norm(x-y))
}

#mixture of truncated normal distribution
rmix_truncnorm <- function(n=1,left_prob=0.5,a=1,b=Inf,mean=0,sd=1)
{
  x <- NULL
  for (i in 1:n)
  {
    tau <- rbinom(n=1,size=1,prob = left_prob) 
    c <- (1-tau)*rtruncnorm(n=1,a=(-1*b),b=(-1*a),mean=(-1*mean),sd=sd) + tau*rtruncnorm(n=1,a=a,b=b,mean=mean,sd=sd)
    x <- append(x,c)
  }
  return(x)
}

#calculate MSE
cal_NMSE <- function(y_te, X_te, w_hat)
{
  n <- length(y_te)
  return(sqL2_norm(y_te-X_te%*%w_hat)/(n-1)/var(y_te))
}

#simulation_studies
# make_dataset <- function(T,n_k,p,C,C_sd,v_var,o_var,kappa=0,ep_sd,rho_s=0,rho_o=0,rowwise_out=FALSE,ex_random=FALSE,low_a=1,up_b=Inf,trunc_mean= 3)
# {
#   
#   h <- T/C 
#   
#   W_c <- matrix(0,C,p) 
#   W_k <- matrix(0,T,p) 
#   O_k <- matrix(0,T,p) 
#   
#   W_ast <- matrix(0,T,p) 
#   
#   clus_split_ind <- sample.int(C, size = p, replace = TRUE , prob = NULL) 
#   
#   for (i in 1:C)
#   {
#     W_c[i,which(clus_split_ind==i)] <- 1
#     for(j in 1:p)
#     {
#       if(W_c[i,j]==1)
#         W_c[i,j] <- rnorm(1,0,sd=C_sd)
#     }
#   }
#   U_ast_list <- vector("list",length = T)
#   for (i in 0:(C-1))
#   {
#     for (l in 1:h)
#     {
#       j <- i*h + l
#       U_ast_list[[j]] <- as.vector(W_c[i+1,]) 
#     }
#   }
#   
# 
#   outlier_vec <- rep(0,T)
#   for (i in 0:(C-1)) 
#   {
#     f <- i*h+1
#     l <- (i+1)*h
#     for (a in f:l)
#     {
#       W_k[a,which(clus_split_ind==(i+1))] <- 1
#       outlier_vec[a] <- rbinom(n=1,size=1,prob = kappa) #確立kappaでa版目のタスクは外れ値
#       if(outlier_vec[a]==1 && rowwise_out==TRUE)
#         O_k[a,] <- mvrnorm(1,rep(0,p),Sigma=Sigma_o)
#       for (j in 1:p) 
#       {
#         if(W_k[a,j]==1)
#         {
#           W_k[a,j] <- rnorm(1,0,sd=sqrt(w_var))
#           if(outlier_vec[a]==1 && rowwise_out==FALSE)
#           {
#             tau <- rbinom(n=1,size=1,prob = 0.5) 
#             O_k[a,j] <-   (1-tau)*rtruncnorm(-1,a=(-1*up_b),b=(-1*low_a),mean=trunc_mean,sd=sqrt(o_var)) + tau*rtruncnorm(1,a=low_a,b=up_b,mean=trunc_mean,sd=sqrt(o_var))
#           }
#         }
#       }
#       W_ast[a,] <- W_k[a,] + W_c[i+1,] + O_k[a,]
#     }
#   }
#   
# 
#   W_ast_list <- asplit(W_ast, MARGIN = 1)
#   y_list <- vector("list", length = T)
#   X_list <- vector("list", length = T)
#   n_list <- vector("list", length = T)
#   sigma_list <- vector("list", length = T)
#   
#   
#   
#   Sigma <- matrix(0,p,p)
#   for (i in 1:p)
#   {
#     for(j in 1:p)
#     {
#       Sigma[i,j] = rho_s^abs(i-j)
#     }
#   }
#   
#   if(ex_random==TRUE)
#   {
#     sigma_vec <- abs(rnorm(T,mean=10,sd=10))
#   }
#   else
#   {
#       #sigma_vec <-sample(seq(1,5,2),T,replace = TRUE)
#     sigma_vec <- rep(ep_sd,T)
#   }
# 
#   for (k in 1:T)
#   {
#     n_list[[k]] <- n_k
#     sigma_list[[k]] <- abs(rnorm(1,mean=10,sd=10))
#     X_list[[k]] <- mvrnorm(n=n_list[[k]],mu=rep(0,p),Sigma = Sigma)
#     y_list[[k]] <- X_list[[k]]%*%W_ast_list[[k]] + rnorm(n_list[[k]],mean = 0,sigma_vec[k])
#   }
#   
#   true_member <- rep(1:C, times = rep(h,C))
#   
#   return(list(X_list,y_list,W_ast_list,U_ast_list,sigma_vec,true_member))
# }


#make_dataset_in_simulation_studies_for_Case1
make_dataset_outlier_case1 <- function(T,n_k,p,C,C_sd,v_var,o_var,kappa=0,ep_sd,rho_s=0,low_a=1,up_b=Inf,trunc_mean= 3,rowwise_out=FALSE)
{
  h <- T/C #number of task in a cluster
  
  W_c <- matrix(0,C,p) 
  W_k <- matrix(0,T,p) 
  O_k <- matrix(0,T,p)
  
  W_ast <- matrix(0,T,p) #真の回帰係数
  
  clus_split_ind <- sample.int(C, size = p, replace = TRUE , prob = NULL) #変数をC個のクラスタに分割
  #クラスタ毎の重心を与える
  for (i in 1:C)
  {
    W_c[i,which(clus_split_ind==i)] <- 1
    for(j in 1:p)
    {
      if(W_c[i,j]==1)
        W_c[i,j] <- rnorm(1,0,sd=C_sd)
    }
  }
  U_ast_list <- vector("list",length = T)
  for (i in 0:(C-1))
  {
    for (l in 1:h)
    {
      j <- i*h + l
      U_ast_list[[j]] <- as.vector(W_c[i+1,])
    }
  }
  
  
  outlier_vec <- rep(0,T)
  for (i in 0:(C-1)) 
  {
    f <- i*h+1
    l <- (i+1)*h
    for (a in f:l)
    {
      W_k[a,which(clus_split_ind==(i+1))] <- 1
      outlier_vec[a] <- rbinom(n=1,size=1,prob = kappa)
      if(outlier_vec[a]==1 && rowwise_out==TRUE)
        O_k[a,] <- rmix_truncnorm(n=p, a=low_a,b=up_b,mean=trunc_mean,sd=sqrt(o_var))
      for (j in 1:p) 
      {
        if(W_k[a,j]==1)
        {
          W_k[a,j] <- rnorm(1,0,sd=sqrt(w_var))
          if(outlier_vec[a]==1 && rowwise_out==FALSE)
          {
            O_k[a,j] <- rmix_truncnorm(n=1, a=low_a,b=up_b,mean=trunc_mean,sd=sqrt(o_var))
          }
        }
      }
      W_ast[a,] <- W_k[a,] + W_c[i+1,] + O_k[a,]
    }
  }
  
  outlier_ind <- which(outlier_vec==1)
  
  W_ast_list <- asplit(W_ast, MARGIN = 1)
  y_list <- vector("list", length = T)
  X_list <- vector("list", length = T)
  n_list <- vector("list", length = T)
  sigma_list <- vector("list", length = T)
  
  
  
  Sigma <- matrix(0,p,p)
  for (i in 1:p)
  {
    for(j in 1:p)
    {
      Sigma[i,j] = rho_s^abs(i-j)
    }
  }
  
  sigma_vec <- rep(ep_sd,T)
  
  for (k in 1:T)
  {
    n_list[[k]] <- n_k
    sigma_list[[k]] <- abs(rnorm(1,mean=10,sd=10))
    X_list[[k]] <- mvrnorm(n=n_list[[k]],mu=rep(0,p),Sigma = Sigma)
    y_list[[k]] <- X_list[[k]]%*%W_ast_list[[k]] + rnorm(n_list[[k]],mean = 0,sigma_vec[k])
  }
  
  true_member <- rep(1:C, times = rep(h,C))
  
  return(list(X_list,y_list,W_ast_list,U_ast_list,sigma_vec,true_member,outlier_ind))
}

#make_dataset_in_simulation_studies_for_Case2
make_dataset_outlier_case2 <- function(T,n_k,p,C,C_sd,v_var,kappa=0,ep_sd,rho_s=0,rowwise_out=FALSE)
{
  
  h <- T/C   
  W_c <- matrix(0,C,p)
  W_k <- matrix(0,T,p)
  O_k <- matrix(0,T,p)
  
  W_ast <- matrix(0,T,p)
  
  clus_split_ind <- sample.int(C, size = p, replace = TRUE , prob = NULL)
  
  for (i in 1:C)
  {
    W_c[i,which(clus_split_ind==i)] <- 1
    for(j in 1:p)
    {
      if(W_c[i,j]==1)
        W_c[i,j] <- rnorm(1,0,sd=C_sd)
    }
  }
  U_ast_list <- vector("list",length = T)
  for (i in 0:(C-1))
  {
    for (l in 1:h)
    {
      j <- i*h + l
      U_ast_list[[j]] <- as.vector(W_c[i+1,]) 
    }
  }
  
  
  outlier_vec <- rep(0,T)
  for (i in 0:(C-1)) 
  {
    f <- i*h+1
    l <- (i+1)*h
    for (a in f:l)
    {
      W_k[a,which(clus_split_ind==(i+1))] <- 1
      outlier_vec[a] <- rbinom(n=1,size=1,prob = kappa) 
      for (j in 1:p) 
      {
        if(W_k[a,j]==1)
        {
          if(outlier_vec[a]==1 && rowwise_out==FALSE)
          {
            W_ast[a,j] <-  runif(1,min=-10,max=10)
          }
          else W_k[a,j] <- rnorm(1,0,sd=sqrt(w_var))
        }
      }
      if(outlier_vec[a]==1&&rowwise_out==TRUE)
        W_ast[a,] <-  runif(p,min=-10,max=10)
      else
      {
        W_ast[a,] <- W_k[a,] + W_c[i+1,]
      }
    }
  }
  
  
  W_ast_list <- asplit(W_ast, MARGIN = 1)
  y_list <- vector("list", length = T)
  X_list <- vector("list", length = T)
  n_list <- vector("list", length = T)
  sigma_list <- vector("list", length = T)
  
  
  
  Sigma <- matrix(0,p,p)
  for (i in 1:p)
  {
    for(j in 1:p)
    {
      Sigma[i,j] = rho_s^abs(i-j)
    }
  }
  
  
  sigma_vec <- rep(ep_sd,T)
  
  for (k in 1:T)
  {
    n_list[[k]] <- n_k
    X_list[[k]] <- mvrnorm(n=n_list[[k]],mu=rep(0,p),Sigma = Sigma)
    y_list[[k]] <- X_list[[k]]%*%W_ast_list[[k]] + rnorm(n_list[[k]],mean = 0,sigma_vec[k])
  }
  
  outlier_ind <- which(outlier_vec==1)
  true_member <- rep(1:C, times = rep(h,C))
  
  return(list(X_list,y_list,W_ast_list,U_ast_list,sigma_vec,true_member,outlier_ind))
}


#output_
make_real_dataset <- function(tr_prop=0.7,v_prop=0.2,te_prop=0.1,data="school",path)
{
  
  if(data=="school")
  {
    T = 139
    y_list <- vector("list",length =T)
    X_list <- vector("list",length =T)
    d_list <- vector("list",length =T)
    
    for (k in 1:T)
    {
      k_char <- as.character(k)
      yfile <- paste(path,"school_Y",k_char,".csv",sep="")
      Xfile <- paste(path,"school_X",k_char,".csv",sep="")
      y <- read.csv(yfile,header = F)
      X <- read.csv(Xfile,header = F)
      d_list[[k]] <- cbind(y,X)
    }
    
    #split_samples_into_samples_for_test_and_learning
    d_split <- lapply(d_list,prop=1-te_prop,initial_split)
    d_list_le <- lapply(d_split, training)
    d_list_te <- lapply(d_split,testing)
    
    #split_samples_for_learning_into_samples_for_training_and_validation
    d_split_le <- lapply(d_list_le,prop=1-v_prop/(tr_prop+v_prop),initial_split)
    d_list_tr <- lapply(d_split_le, training)
    d_list_v <- lapply(d_split_le, testing)
    
    #make_response_and_explanatory_variables
    y_list_tr <- d_list_tr %>% purrr::map(~.x[,1])
    X_list_tr <- d_list_tr %>% purrr::map(~.x[,-1])
    
    y_list_v <- d_list_v %>% purrr::map(~.x[,1])
    X_list_v <- d_list_v %>% purrr::map(~.x[,-1])
    
    y_list_te <- d_list_te %>% purrr::map(~.x[,1])
    X_list_te <- d_list_te %>% purrr::map(~.x[,-1])
    
    
    y_list_tr <-lapply(y_list_tr,as.vector)
    y_list_v <-lapply(y_list_v,as.vector)
    y_list_te <-lapply(y_list_te,as.vector)
    
    

    X_list_tr <-lapply(X_list_tr,as.matrix) 
    X_list_v <-lapply(X_list_v,as.matrix) 
    X_list_te <-lapply(X_list_te,as.matrix)
    
    #
    tr_data <-list(y_list_tr,X_list_tr)
    v_data <-list(y_list_v,X_list_v)
    te_data <-list(y_list_te,X_list_te)
  }
  
  if(data=="landmine")
  {
    T = 29
    y_list <- vector("list",length =T)
    X_list <- vector("list",length =T)
    d_list <- vector("list",length =T)
    down_d_list <- vector("list",length =T)
    
    
    for (k in 1:T)
    {
      k_char <- as.character(k)
      yfile <- paste(path,"landmine_Y",k_char,".csv",sep="")
      Xfile <- paste(path,"landmine_X",k_char,".csv",sep="")
      y <- read.csv(yfile,header = F,encoding='Shift_JIS')
      X <- read.csv(Xfile,header = F,encoding='Shift_JIS')
      # X <- read.csv(Xfile,header = F,encoding='Shift_JIS')
      d_list[[k]] <- cbind(y,X)
      colnames(d_list[[k]])[1] <- "Res"
      #down_sampling
      down_d_list[[k]] <- rbind(d_list[[k]][sample(which(d_list[[k]][,1]==0),sum(d_list[[k]][,1])),],d_list[[k]][which(d_list[[k]][,1]==1),])
    }
    
    #lapply(down_d_list, scale , center=TRUE, scale=TRUE)
    
    
    #split_samples_into_samples_for_test_and_learning
    d_split <- lapply(down_d_list,prop=1-te_prop,strata = "Res",initial_split)
    d_list_le <- lapply(d_split, training)
    d_list_te <- lapply(d_split,testing)
    
    #split_samples_for_learning_into_samples_for_training_and_validation
    d_split_le <- lapply(d_list_le,prop=1-v_prop/(tr_prop+v_prop),strata = "Res",initial_split)
    d_list_tr <- lapply(d_split_le, training)
    d_list_v <- lapply(d_split_le, testing)
    
    #make_response_and_explanatory_variables
    y_list_tr <- d_list_tr %>% purrr::map(~.x[,1])
    X_list_tr <- d_list_tr %>% purrr::map(~.x[,-1])
    
    y_list_v <- d_list_v %>% purrr::map(~.x[,1])
    X_list_v <- d_list_v %>% purrr::map(~.x[,-1])
    
    y_list_te <- d_list_te %>% purrr::map(~.x[,1])
    X_list_te <- d_list_te %>% purrr::map(~.x[,-1])
    
    
    y_list_tr <-lapply(y_list_tr,as.vector)
    y_list_v <-lapply(y_list_v,as.vector)
    y_list_te <-lapply(y_list_te,as.vector)
    
    
    X_list_tr <-lapply(X_list_tr,data.matrix)
    X_list_v <-lapply(X_list_v,data.matrix)
    X_list_te <-lapply(X_list_te,data.matrix)
    
    
    #standardize_X
    X_list_tr <-lapply(X_list_tr,original_scale) 
    X_list_v <-lapply(X_list_v,original_scale)
    X_list_te <-lapply(X_list_te,original_scale)
    
    
    tr_data <-list(y_list_tr,X_list_tr)
    v_data <-list(y_list_v,X_list_v)
    te_data <-list(y_list_te,X_list_te)
  }
  data_list <- list(tr_data,v_data,te_data)
  return(data_list)
}
