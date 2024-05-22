
prox <- function(Z,lambda_vec)
{
  norm_vec <- apply(Z, MARGIN = 1, FUN = L2_norm)
  norm_vec[norm_vec<1e-10] <- 1e-5 #to avoid error
  return(pmin(lambda_vec,norm_vec)*Z/norm_vec)
}

#thresholding_function_grouplasso
g.soft.th = function(a,lambda,gamma=0)
{
  return(max(0,1-lambda/L2_norm(a))*a)
}

#thresholding_function_groupMCP
gMCP = function(a, lambda,gamma)
{
  if (L2_norm(a)<=gamma*lambda)
  {
    return(gamma/(gamma-1)*g.soft.th(a,lambda))
  }
  else return(a)
}

#thresholding_function_groupSCAD
gSCAD = function(a,lambda,gamma)
{
  b <- L2_norm(a)
  if(b<=2*lambda)
  {
    return(g.soft.th(a,lambda))
  }
  else
  {
    if(b <= gamma*lambda)
    {
      return((gamma-1)/(gamma-2)*g.soft.th(a,gamma*lambda/(gamma-1)))
    }
    else
    {
      return(a)
    }
  }
}

#thresholding_function_groupHard_ridge
Hard_ridge = function(a,lambda,gamma)
{
  p <- length(a)
  if(L2_norm(a)<lambda)
  {
    return(rep(0,p))
  }
  else return(a/(1+gamma))
}

#transforming_Rmatrix_to_vector
RtoVec <- function(R)
{
  R_vec <- 0
  T <- ncol(R)
  for (i in 1:T) 
  {
    R_vec <- append(R_vec,R[i,-(1:i)])
  }
  R_vec <- R_vec[R_vec!=0] #extract_nonzero_elements
  return(R_vec)
}

#transforming_R_to_Aepsilon
make_A <- function(R,R_vec)
{
  T <- ncol(R)
  A_ep <- matrix(0,length(R_vec),T)
  c<-1
  for (i in 1:(T-1)) 
  {
    l <- i+1 
    for (j in l:T)
    {
      if(R[i,j]!=0)
      {
        e <- rep(0,T)
        e[i] <- 1
        e[j] <- -1
        A_ep[c,] <- e
        c <- c+1
      }
    }
  }
  return(A_ep)
}

#calculating_linear_ridge_estimator_with_center_u
centered_ridge <- function(y,solve,Xy,u,lambda1)
{
  n <- length(y)
  u <- as.vector(u)
  return(solve %*% (Xy/n+lambda1*u))
}

#calculating_logistic_ridge_estimator_with_center_u
centered_ridge_logistic <- function(y,X,u,lambda_1,lambda_3 = 0.01,lambda_4=0)
{
  #lambda3 ridge for intercept
  #lambda4 ridge for coefficients of features 
  n <- nrow(X)
  p <- ncol(X)
  X <- cbind(rep(1,n),X)

  I_n <- diag(n)
  I_p1 <- diag(p+1)
  
  Lambda <- lambda_1*I_p1
  Lambda[1,1] <- 0
  Gamma <- lambda_4*I_p1
  Gamma[1,1] <- lambda_3
  pi_vec <- rep(0,n)
  # w <- rep(0,p)
  w <-as.vector(c(0,u))
  u_dash <- c(0,u)
  
  eps <-10000
  Pi <- matrix(0,n,n)
  count <- 0
  while(eps>0.001 && count < 20)
  {
    count <- count + 1
    w_old <- w
    pi_vec_old <- pi_vec
    
    for (i in 1:n)
    {
      pi_vec[i] <- pi(X[i,],w)
    }
    
    
    H <- diag(pi_vec*(1-pi_vec))
    
    w <- as.vector(w + solve(t(X)%*%H%*%X/n+Lambda+Gamma)%*%(as.vector(t(X)%*%(y-pi_vec)/n-(Lambda%*%(w-u_dash)+Gamma%*%w))))
    
    eps <- max(abs(w-w_old))
  }
  return(w)
}



MTLCVX_estimator<- function(y_list,X_list,R,lambda1,lambda2,rho=1,thre_ADMM=0.01,max_iter=200, response="gaussian")
{
  if(response=="gaussian")
    return(MTLCVX_linear(y_list,X_list,R,lambda1,lambda2,rho,thre_ADMM,max_iter))
  if(response=="binomial")
    return(MTLCVX_logistic(y_list,X_list,R,lambda1,lambda2,rho,thre_ADMM,max_iter))
}


MTLRRC_estimator<- function(y_list,X_list,R,lambda1,lambda2,lambda3,gamma=3,thre_type =g.soft.th,rho=1,thre_ADMM=0.01,max_iter=200,response="gaussian")
{
  if(response=="gaussian")
    return(MTLRCC_linear(y_list,X_list,R,lambda1,lambda2,lambda3,gamma,thre_type=thre_type,  rho,thre_ADMM,max_iter))
  if(response=="binomial")
    return(MTLRCC_logistic(y_list,X_list,R,lambda1,lambda2,lambda3,gamma,thre_type =thre_type, rho,thre_ADMM,max_iter))
}


#MTLK
MTLK_estimator <- function(y_list,X_list,W_ini,lambda,k,thre_ADMM=0.1,max_iter=200,response="gaussian")
{
  if(response=="gaussian")
    return(MTLK_linear(y_list,X_list,W_ini,lambda,k,thre_ADMM,max_iter))
  if(response=="binomial")
    return(MTLK_logistic(y_list,X_list,W_ini,lambda,k,thre_ADMM,max_iter))
}

MTLCVX_linear <- function(y_list,X_list,R,lambda_1,lambda_2,rho=1,thre_ADMM=0.01,max_iter=200)
{
  T <- length(y_list)
  p <- ncol(X_list[[1]])
  #regression coefficients and centroids
  W <- matrix(0,T,p)
  U <- matrix(0,T,p)
  I_p <- diag(p)
  R_vec <- RtoVec(R)
  A_ep  <- make_A(R,R_vec)
  G <- t(A_ep) %*% A_ep
  #stepsize in FISTA
  eta <- 1/(lambda_1+2*rho*max(diag(G)))
  #prior calculation
  solve_list <- vector("list",length = T)
  for (k in 1:T)
  {
    solve_list[[k]] <- solve(t(X_list[[k]])%*%X_list[[k]]/nrow(X_list[[k]])+lambda_1*I_p)
  }
  Xy_list <- mapply(norm,X_list,y_list,SIMPLIFY = FALSE)
  
  #Lagrangian multiplier
  S <- matrix(0,length(R_vec),p)
  
  eps <- 10000
  count <- 0
  while(eps>thre_ADMM && count<max_iter)
  {
    count <- count + 1
    old_W <- W
    
    #convert_matrix_to_U
    U_list <- asplit(U,MARGIN = 1)
    #update_W
    W <- t(mapply(centered_ridge, y_list,solve_list, Xy_list,U_list, lambda_1,SIMPLIFY =T)) #回帰係数の更新 通常版
    #update_U
    #initialize FISTA
    alpha_t <- 1
    Z <- matrix(0,T,p)
    eps_F <-10000 
    while (eps_F>0.001)
    {
      Z_old <- Z
      alpha_t_old <- alpha_t
      Z <- U - eta * (lambda_1 * (U-W) + t(A_ep) %*% prox(S + rho * A_ep %*% U,lambda_2 * R_vec))
      alpha_t <- (1+sqrt(1+4*alpha_t_old^2))/2
      U <- Z + (alpha_t_old-1)/alpha_t*(Z-Z_old)
      eps_F <- max(Z-Z_old)
    }
    U <- Z
    #update_S
    S <- prox(S + rho * A_ep %*% U,lambda_2 * R_vec)
    
    #velue for evaluation of convergence
    eps <- max(abs(W-old_W))
  }
  W_list <- asplit(W,MARGIN = 1)
  U_list <- asplit(U,MARGIN = 1)
  output <- list(W_list,U_list)
  return(output)
}

MTLCVX_logistic <- function(y_list,X_list,R,lambda_1,lambda_2,rho=1,thre_ADMM=0.01,max_iter=200)
{
  T <- length(y_list)
  p <- ncol(X_list[[1]])
  #regression coefficients and centroids
  W <- matrix(0,T,p+1)
  U <- matrix(0,T,p)
  
  I_p <- diag(p)
  R_vec <- RtoVec(R)
  A_ep  <- make_A(R,R_vec)
  G <- t(A_ep) %*% A_ep
  #  #stepsize in FISTA
  eta <- 1/(lambda_1+2*rho*max(diag(G)))
  
  #Lagrangian multiplier
  S <- matrix(0,length(R_vec),p)
  
  
  eps <- 10000
  count <- 0
  while(eps>thre_ADMM && count<max_iter)
  {
    count <- count + 1
    old_W <- W
    
    #    #convert_matrix_to_list
    U_list <- asplit(U,MARGIN = 1)
    #update_W
    W <- t(mapply(centered_ridge_logistic, y_list, X_list,U_list, lambda_1,SIMPLIFY =T)) #回帰係数の更新 通常版
    W_outint <- W[,-1]
    #initialize FISTA
    alpha_t <- 1
    Z <- matrix(0,T,p)
    eps_F <-10000 #FISTAの収束
    while (eps_F>0.001)
    {
      Z_old <- Z
      alpha_t_old <- alpha_t
      Z <- U - eta * (lambda_1 * (U-W_outint) + t(A_ep) %*% prox(S + rho * A_ep %*% U,lambda_2 * R_vec))
      alpha_t <- (1+sqrt(1+4*alpha_t_old^2))/2
      U <- Z + (alpha_t_old-1)/alpha_t*(Z-Z_old)
      eps_F <- max(Z-Z_old)
    }
    U <- Z
    #update_S
    S <- prox(S + rho * A_ep %*% U,lambda_2 * R_vec)
    
    #velue for evaluation of convergence
    eps <- max(abs(W-old_W))
  }
  
  W_list <- asplit(W,MARGIN = 1)
  U_list <- asplit(U,MARGIN = 1)
  output <- list(W_list,U_list)
  return(output)
}

MTLRCC_linear <- function(y_list,X_list,R,lambda1,lambda2,lambda3,gamma,thre_type = g.soft.th,rho=1,thre_ADMM=0.01,max_iter=200)
{
  T <- length(y_list)
  p <- ncol(X_list[[1]])
  #regression coefficients, centroids and outlier parameters
  W <- matrix(0,T,p)
  U <- matrix(0,T,p)
  O <- matrix(0,T,p)
  
  Q <- matrix(0,T,p)
  I_p <- diag(p)
  R_vec <- RtoVec(R)
  A_ep  <- make_A(R,R_vec)
  G <- t(A_ep) %*% A_ep
  #stepsize in FISTA
  eta <- 1/(lambda1+2*rho*max(diag(G)))
  
  #prior calculation
  XTX_list <- lapply(X_list,sqL2_norm)
  Xy_list <- mapply(norm,X_list,y_list,SIMPLIFY = FALSE)
  solve_list <- vector("list",length = T)
  for (k in 1:T)
  {
    solve_list[[k]] <- solve(t(X_list[[k]])%*%X_list[[k]]/nrow(X_list[[k]])+lambda1*I_p)
  }
  lam3_lam1 <- lambda3/lambda1
  
  #Lagrangian multiplier
  S <- matrix(0,length(R_vec),p)
  
  eps <- 10000
  count <- 0
  while(eps>0.01 && count<200 )
  {
    count <- count + 1
    old_W  <- W
    
    #convert_matrix_to_list
    U_list <- asplit(U,MARGIN = 1)
    Q_list <- asplit(Q,MARGIN = 1)

    #update_W
    W <- t(mapply(centered_ridge, y_list,solve_list, Xy_list,Q_list, lambda1,SIMPLIFY =T)) 
    
    #update_U
    #initialize FISTA
    alpha_t <- 1
    Z <- matrix(0,T,p)
    eps_F <-10000
    while (eps_F>0.001)
    {
      Z_old <- Z
      alpha_t_old <- alpha_t
      Z <- U - eta * (lambda1 * (U+O-W) + t(A_ep) %*% prox(S + rho * A_ep %*% U,lambda2 * R_vec))
      alpha_t <- (1+sqrt(1+4*alpha_t_old^2))/2
      U <- Z + (alpha_t_old-1)/alpha_t*(Z-Z_old)
      eps_F <- max(Z-Z_old)
    }
    U <- Z
    #update_O
    O <- t(mapply(thre_type,asplit(W-U,MARGIN = 1),lam3_lam1,gamma,SIMPLIFY =T))
    
    Q  <- O + U
    
    #update_S
    S <- prox(S + rho * A_ep %*% U,lambda2 * R_vec)
    eps <- max(abs(W-old_W))
  }
  W_list <- asplit(W,MARGIN = 1)
  U_list <- asplit(U,MARGIN = 1)
  O_list <- asplit(O,MARGIN = 1)
  
  output <- list(W_list,U_list,O_list)
  return(output)
}

MTLRCC_logistic <- function(y_list,X_list,R,lambda1,lambda2,lambda3,gamma,thre_type = g.soft.th,rho=1,thre_ADMM=0.01,max_iter=200)
{
  
  T <- length(y_list)
  p <- ncol(X_list[[1]])
  #regression coefficients, centroids and outlier parameters
  W <- matrix(0,T,p+1)
  U <- matrix(0,T,p)
  O <- matrix(0,T,p)

  Q <- matrix(0,T,p)
  I_p <- diag(p)
  R_vec <- RtoVec(R)
  A_ep  <- make_A(R,R_vec)
  G <- t(A_ep) %*% A_ep
  
  #stepsize in FISTA
  eta <- 1/(lambda1+2*rho*max(diag(G)))
  
  
  
  lam3_lam1 <- lambda3/lambda1
  
  #Lagrangian multiplier
  S <- matrix(0,length(R_vec),p)
  
  eps <- 10000
  count <- 0
  while(eps>thre_ADMM && count<max_iter )
  {
    count <- count + 1
    old_W  <- W
    
    #convert_matrix_to_list
    U_list <- asplit(U,MARGIN = 1)
    Q_list <- asplit(Q,MARGIN = 1)
    
    #update_W
    W <- t(mapply(centered_ridge_logistic, y_list,X_list, Q_list, lambda1,SIMPLIFY =T))
    W_outint <- W[,-1]
    
    #initialize FISTA
    alpha_t <- 1
    Z <- matrix(0,T,p)
    eps_F <-10000 
    while (eps_F>0.001)
    {
      Z_old <- Z
      alpha_t_old <- alpha_t
      Z <- U - eta * (lambda1 * (U+O-W_outint) + t(A_ep) %*% prox(S + rho * A_ep %*% U,lambda2 * R_vec))
      alpha_t <- (1+sqrt(1+4*alpha_t_old^2))/2
      U <- Z + (alpha_t_old-1)/alpha_t*(Z-Z_old)
      eps_F <- max(Z-Z_old)
    }
    U <- Z
    
    #update_O
    O <- t(mapply(thre_type,asplit(W_outint-U,MARGIN = 1),lam3_lam1,gamma,SIMPLIFY =T))
    
    Q  <- O + U
    #update_S
    S <- prox(S + rho * A_ep %*% U,lambda2 * R_vec)
    eps <- max(abs(W-old_W))
  }
  W_list <- asplit(W,MARGIN = 1)
  U_list <- asplit(U,MARGIN = 1)
  O_list <- asplit(O,MARGIN = 1)
  
  output <- list(W_list,U_list,O_list)
  return(output)
}

MTLk_linear <- function(y_list,X_list,W_ini,lambda,k,thre_ADMM=0.01,max_iter=20)
{
  T <- length(y_list)
  p <- ncol(X_list[[1]])
  #regression coefficients and cluster center
  W <- W_ini
  U <- matrix(0,T,p)
  I_p <- diag(p)

  #prior calculation
  solve_list <- vector("list",length = T)
  for (m in 1:T)
  {
    solve_list[[m]] <- solve(t(X_list[[m]])%*%X_list[[m]]/nrow(X_list[[m]])+lambda*I_p)
  }
  Xy_list <- mapply(norm,X_list,y_list,SIMPLIFY = FALSE)
  
  eps <- 10000
  count <- 0
  while(eps>thre_ADMM && count<max_iter)
  {
    count <- count + 1
    old_W <- W
    #update_U
    km <- kmeans(W,centers=k)
    for (m in 1:T)
    {
      U[m,] <- km$center[km$cluster[m],]
    }
    
    #convert_matrix_to_list
    U_list <- asplit(U,MARGIN = 1)
    
    #update_W
    W <- t(mapply(centered_ridge, y_list,solve_list, Xy_list,U_list, lambda,SIMPLIFY =T)) 
 
    #収束の確認
    eps <- max(abs(W-old_W))
  }
  
  W_list <- asplit(W,MARGIN = 1)
  U_list <- asplit(U,MARGIN = 1)
  
  output <- list(W_list,U_list)
  return(output)
}

MTLk_logistic <- function(y_list,X_list,W_ini,lambda,k,thre_ADMM=0.01,max_iter=200)
{
  T <- length(y_list)
  p <- ncol(X_list[[1]])
  
  #regression coefficients and cluster center
  W <- W_ini
  U <- matrix(0,T,p)
  I_p <- diag(p)
  
  eps <- 10000
  count <- 0
  while(eps>thre_ADMM && count<max_iter)
  {
    count <- count + 1
    old_W <- W
    #update_U
    km <- kmeans(W[,-1],centers=k)
    for (m in 1:T)
    {
      U[m,] <- km$center[km$cluster[m],]
    }
    #convert_matrix_to_list
    U_list <- asplit(U,MARGIN = 1)
    #update_W
    W <- t(mapply(centered_ridge_logistic, y_list,X_list, U_list, lambda,SIMPLIFY =T))
    eps <- max(abs(W-old_W))
  }
  W_list <- asplit(W,MARGIN = 1)
  U_list <- asplit(U,MARGIN = 1)
  
  output <- list(W_list,U_list)
  return(output)
}




