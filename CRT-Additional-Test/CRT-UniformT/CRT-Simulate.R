
source("../Gof/GoF-Settings.R")

## Function to simulate x
simulate_data <- function(Sigma, n) {
  X <- mvrnorm(n, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
  return(X)
}

simulate_data_t <- function(n,p,coefs=c(0.2,0.2,-0.2,-0.2,0.2,0.2),df=3, noise_scale=0.7) {
  
  k=length(coefs)
  X <- matrix(0, n, p)
  X[, 1] <- rt(n, df)  # First column from t-distribution
    
    for (j in 2:p) {
      kj <- min(j - 1, k)  # Number of past columns to use
      beta <- coefs[1:kj]
      prev_cols <- 
      
      # Generate centered linear combination
      lin_comb <- X[, (j - kj):(j - 1), drop = FALSE] %*% beta[1:kj]
      noise <- rt(n, df) * noise_scale
      X[, j] <- lin_comb + noise
    }
    
    return(X)
  }
  

## Function to simulate y
random_beta =function(p){ 
  set.seed(2023)
  beta=rep(0,p)
  p_tilde = min(p,20)
  beta[1:p_tilde] <- runif(p_tilde, 1,2)/sqrt(p_tilde)
  
  return(beta)
}

bq3=function(vals, lq=qnorm(1/3), uq=qnorm(2/3), w=1) {
  res <- vals
  res[vals > uq] <- w
  res[vals >= lq & vals <= uq] <- -2 * w
  res[vals < lq] <- w
  return(res)
}

non_linear_link_regression=function(x, theta){
  basis=cbind(x[,1]^2/2,
              1/(1+x[,2]^2),
              cos(pi*x[,3]),
              sin(pi*x[,4]),
              x[,5]*x[,6], 
              sin( x[,7]*x[,8] ), 
              sin(pi*x[,9])/(4+x[,10]^2), 
              
              x[,11]^2/2,
              1/(1+x[,12]^2),
              cos(pi*x[,13]),
              sin(pi*x[,14]),
              x[,15]*x[,16], 
              sin( x[,17]*x[,18] ), 
              sin(pi*x[,19])/(4+x[,20]^2)
  ) 
  xb = basis %*% c(rep(theta,7),rep(1,7)) 
  return(c(xb))
}


non_linear_link_classification=function(x, theta){
  basis=apply(x[,1:20],2,bq3) / 4
  xb = basis %*% c(rep(theta,10),rep(1,10)) 
  return(c(xb))
}


simulate_y <- function(x,setting,setT,theta,beta,
                       noise_unif=NULL){
  n=nrow(x)
  p=ncol(x)
  if(is.null(noise_unif)){
    noise_unif=runif(n)
  }
  noise_norm=qnorm(noise_unif)
  
  setting_model=substr(setting,start = 3,stop=6)
  if (setting_model == "l_w") {
    scale_beta=beta;
    scale_beta[setT]=scale_beta[setT]*theta
    
    xb = x%*% scale_beta
    y <- xb +  noise_norm
  }  else if (setting_model == "gl_w") {
    scale_beta=beta;
    scale_beta[setT]=scale_beta[setT]*theta
    
    xb = x%*% scale_beta
    mu = binomial()$linkinv(xb)
    y = 1 * (noise_unif<mu) 
  } else if (setting_model == "l_m") {
    xb=non_linear_link_regression(x,theta)
    y= xb +  noise_norm
    
  }else if(setting_model == "gl_m"){
    
    xb=non_linear_link_classification(x, theta)
    mu = binomial()$linkinv(xb)
    y = 1 * (noise_unif<mu)
    
  }else {
    stop("Invalid setting value.")  
  }
  return(y)
}