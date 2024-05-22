Ini_normal <-function(y,X)
{
  return(solve(t(X)%*%X)%*%t(X)%*%y)
}


Ini_glmnet <- function(y,X,response="gaussian",alpha=1)
{
  y <- as.vector(y)
  X <- as.matrix(X)
  p <- ncol(X)

  if(response=="binomial")
  {
    lasso.model.cv <- cv.glmnet(x =X , y =y , family = response, alpha = alpha)
    lasso.model <- glmnet(x = X, y = y, family = response, lambda =lasso.model.cv$lambda.min , alpha = alpha)
    return(as.vector(coef(lasso.model)))
  }
  
  if(response=="gaussian")
  {
    lasso.model.cv <- cv.glmnet(x =X , y =y , family = response, alpha = alpha)
    lasso.model <- glmnet(x = X, y = y, family = response, lambda =lasso.model.cv$lambda.min , alpha = alpha)
    
    return(as.vector(lasso.model$beta))
  }
}

R_kNN <- function(k,X)
{
  dis_mat <- distance_matrix(X)
  n <- nrow(dis_mat)
  S <- matrix(0,n,n)
  for(i in 1:n)
  {
    dis_mat[i,i] <- 10000000
  }
  for(i in 1:n)
  {
    for(j in 1:k)
    {
      ind <- which(min(dis_mat[i,])==dis_mat[i,])
      S[i,ind] <- 1
      dis_mat[i,ind] <- 100000
    }
  }
  
  
  
  R <- (S + t(S))/2
  return(R)
}

R_gausskernel<-function(k,phi,X)
{
  
  dis_mat <- distance_matrix(X)
  L2_dis <- dis_mat
  n <- nrow(dis_mat)
  T <- matrix(0,n,n)
  R <- matrix(0,n,n)
  for(i in 1:n)
  {
    dis_mat[i,i] <- 10000000
  }
  for(i in 1:n)
  {
    for(j in 1:k)
    {
      ind <- which(min(dis_mat[i,])==dis_mat[i,])
      T[i,ind] <- 1
      dis_mat[i,ind] <- 100000
    }
  }
  for (i in 1:n)
  {
    for(i_dash in 1:n)
    {
      R[i,i_dash] <- T[i,i_dash]*exp(-1*phi*L2_dis[i,i_dash]^2)
    }
  }
  return((t(R)+R)/2)
}

R_unif= function(n)
{
  R <- matrix(1,n,n)
  return(R-diag(n))
}

distance_matrix_scaled <- function(X)
{
  n <- nrow(X)
  scale_X <- scale(X)
  d_X <- matrix(0,n,n)
  
  for(i in 1:n)
  {
    for(i_dash in 1:n)
    {
      d_X[i,i_dash] <-sqrt(vec_norm(scale_X[i,]-scale_X[i_dash,]))
    }
  }
  return(d_X)
}

distance_matrix <- function(X)
{
  n <- nrow(X)
  d_X <- matrix(0,n,n)
  
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      d_X[i,j] <- L2_norm(X[i,]-X[j,])
    }
  }
  return(d_X)
}

#true each element of R is obtained in probability P_{R} (=prob)
R_prob <- function(T,C,prob=1.0)
{
  R <- matrix(0,T,T)
  c <- 0
  h <- T/C
  for (kth_clus in 1:C)
  {
    clus_f_ind <- (kth_clus-1)*h+1
    clus_l_ind <- kth_clus*h
    for (i in clus_f_ind:clus_l_ind)
    {
      g <- i + 1
      if(g>T)
      {
        return(R)
      }
      for (j in g:T)
      {
        if(j<=clus_l_ind)
        {
          R[i,j] <- sample(x=c(1,0),size=1,prob=c(prob,1-prob))
          R[j,i] <- R[i,j] 
        }
        else
        {
            R[i,j] <- sample(x=c(0,1),size=1,prob=c(prob,1-prob))
            R[j,i] <- R[i,j]
        }
      }
    }
  }
  # for (i in 1:T)
  # {
  #   R[i,i] <- 0
  # }
  
  return(R)
}

