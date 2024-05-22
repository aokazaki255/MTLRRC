norm <- function(x,y)
{
  return(t(x)%*%y)
}

#L2ノルムの二乗を計算
sqL2_norm <- function(x)
{
  return(t(x)%*%x)
}

#L2ノルムを計算
L2_norm <- function(x)
{
  return(sqrt(t(x)%*%x))
}

#フロベニウスノルムの計算
F_norm <- function(X)
{
  return(sqrt(sum(diag(X%*%t(X)))))
}

vec_sum <- function(X)
{
  sumX <- apply(X,2,sum)
  return(sumX)
}

hotelingT2 <- function(X,q=0.95,dis="xi")
{
  #ホテリングのT2法
  n <- nrow(X)
  p <- ncol(X)

  X_centered <- scale(X,center = TRUE,scale = FALSE)
  omega <- var(X)
  omega_inv <- solve(omega)
  
  hoteling_t <- rep(0,n)
  
  for (i in 1:n)
  {
    hoteling_t[i] <- t(X_centered[i,])%*%omega_inv%*%X_centered[i,]
  }
  
  #x^2の95%点
  if(dis=="xi")
  {
    th_95 <- qchisq(q, df=p, ncp=0, lower.tail = TRUE, log.p = FALSE)
  }
  if(dis=="F")
  {
    th_95 <- (n-1)*p/(n*(n-p)) * qf(q, df1=p,df2 =(n-p),ncp=0, lower.tail = TRUE, log.p = FALSE)
  }
  #F分布の95%点
  est_outlier_ind <- which(hoteling_t>=th_95)
  return(est_outlier_ind)
} 

cal_Vlikelihood <- function(y_v,X_v,w)
{
  LH <- 0
  n <- length(y_v)
  
  y_hat_v <- rep(0,n)
  lh_v <- rep(0,n)
  for (i in 1:n)
    {
      y_hat_v[i] <- pi(c(1,X_v[i,]),w)
      lh_v[i] <- log((y_hat_v[i]^y_v[i])*((1-y_hat_v[i])^(1-y_v[i])))
    }
  return(-1*mean(lh_v))
}

pi <- function(x,w)
{
  return(1-1/(1+exp(t(x)%*%w)))
}

original_scale <-function(X)
{
 means <-  apply(X, MARGIN =2,FUN="mean")
 sds <-  apply(X, MARGIN =2,FUN="sd")
 centered  <- sweep(X,2,means)
 return(sweep(centered,2,sds,FUN="/"))
}