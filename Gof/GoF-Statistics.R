# GoF STATISTICS --------------------------------------------------------------------
library(Rfast)


## Residual correlation based methods -------

# Helper function to get residuals
# Assume x has been de-mean
get_residuals <- function(x,i, neighbors) {
  if (length(neighbors) == 0) {
    return(x[,i])
  } else {
    model = lmfit(x[, neighbors], x[,i])
    return(model$residuals)
  }
}

# Assume intercept is needed
get_residuals_intercept <- function(x,i, neighbors) {
  if (length(neighbors) == 0) {
    return(x[,i] - mean(x[,i]))
  } else {
    model = lmfit(cbind(1,x[, neighbors]), x[,i])
    return(model$residuals)
  }
}

## 
pvalue_stablizing_cor=function(rho,df){
  z_stat = sqrt(df) * log( (1+rho) / (1 - rho)) / 2
  pv = 2 * pnorm(abs(z_stat), lower.tail=FALSE)
  return(pv)
}

# Compute the Pairwise Residual Correlation (PRC) p-values
PRC_pvalues <- function(x, graph){
  
  p = ncol(x)
  n = nrow(x)
  
  pv.exact = matrix(1, nrow = p, ncol = p)
  pv.stab = matrix(1, nrow = p, ncol = p)
  
  
  if( 2*max(colSums(graph)) + 3 > n){
    ## print('Null graph too large for PRC!')
  }else{
    
    
    x=scale(x,center=T,scale=F)   # de-mean to avoid cbind when getting residual
    
    # Find neighbors for each variable based on the adjacency graph
    neig_list = lapply(1:p, function(i) setdiff(which(graph[i, ] != 0), i))
    
    # Loop over all pairs of variables
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        if (graph[i, j] == 0){  # if i and j are not adjacent
          
          # Find the union of neighbors for i and j
          neig = union(neig_list[[i]], neig_list[[j]])
          # Calculate residuals
          res_i = get_residuals(x, i, neig)
          res_j = get_residuals(x, j, neig)
          
          # Compute residual correlation and test statistics
          res_cor = cor(res_i,res_j)
          df = n - length(neig) - 2
          
          # Compute p-value 
          t_stat = sqrt(df) * res_cor / sqrt(1 - res_cor^2)
          pv.exact[i, j] = 2 * pt(abs(t_stat), df, lower.tail=FALSE)
          pv.stab[i,j]=pvalue_stablizing_cor(res_cor,df)
          
          
          
        }else {
          pv.exact[i, j] = pv.stab[i,j] = 1
        }
        # Store symmetric values
        pv.exact[j, i] = pv.exact[i, j]
        pv.stab[j, i] = pv.stab[i, j] 
        
      }
    }
  }
  
  return(list(pv.exact=pv.exact,pv.stab=pv.stab))   
}


# Compute the Efficient Residual Correlation (ERC) p-values
ERC_pvalues <- function(x, graph){
  p = ncol(x)
  n = nrow(x)
  pv.exact = matrix(1, nrow = p, ncol = p)
  pv.stab = matrix(1, nrow = p, ncol = p)
  
  
  
  if( max(colSums(graph)) + 3 > n){
    ## print('Null graph too large for ERC!')
  }else{
    
    x=scale(x,center=T,scale=F)   # de-mean to avoid cbind when getting residual
    
    
    # Create a list to store neighbors for each variable
    residual = matrix(NA,nr=n,nc=p)
    size_neig=c()
    
    # Do linear regression on neighbors for each variable based on the adjacency graph
    for (i in 1:p){
      i_neig = setdiff(which(graph[i,] != 0), i)
      residual[,i]=get_residuals(x, i,i_neig)
      size_neig[i]=length(i_neig)
    }
    
    # Loop over all pairs of variables
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        if (graph[i, j] == 0){  # if i and j are not adjacent
          
          # Compute residual correlation and test statistics
          res_i = residual[,i]
          res_j = residual[,j]
          
          res_cor = cor(res_i,res_j)
          df=n-2-min(size_neig[i],size_neig[j])
          
          pv.exact[i, j] = pbeta(abs(res_cor),shape1=1,shape2=df,lower.tail = F)
          pv.stab[i, j] = pvalue_stablizing_cor(res_cor,df)
          
          
        } else {
          pv.exact[i, j] = pv.stab[i,j] = 1
        }
        # Store symmetric values
        pv.exact[j, i] = pv.exact[i, j]
        pv.stab[j, i] = pv.stab[i, j] 
      }
    }
  }
  return(list(pv.exact=pv.exact,pv.stab=pv.stab))   
}

## compute one statistic
Stat_RC=function(x, graph,
                 residual_type="PR",RefDist="exact",
                 moment_order=1,
                 threshold=0.1){
  p=ncol(x)
  
  if(residual_type=='PR'){
    pv_list=PRC_pvalues(x, graph)
  }else if(residual_type=='ER'){
    pv_list=ERC_pvalues(x,graph)
  }else {
    stop("Unknown residual type.")
  }
  
  if(RefDist=='exact'){
    pv=pv_list$pv.exact
  }else if (RefDist=='stablizing') {
    pv=pv_list$pv.stab
  }else{
    stop("Unknown reference distribution.")
  }
  
  
  
  z_score = - qnorm(pv/2)
  
  # Calculate the final statistics based on the threshold
  s = 0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if ( graph[i,j]==0 & pv[i, j] < threshold){
        s = s + abs(z_score[i, j])^moment_order
      }
    }
  }
  return(s)
}

# compute all statistics at once for simulation studies
All_Stat_RC=function(x,graph,residual_type="PR",
                     threshold=0.1){
  p=ncol(x)
  if(residual_type=='PR'){
    pv_list=PRC_pvalues(x, graph)
  }else if(residual_type=='ER'){
    pv_list=ERC_pvalues(x,graph)
  }else {
    stop("Unknown residual type.")
  }
  
  all_rc_stat=c()
  
  for(i in 1:2){
    pv=pv_list[[i]] 
    
    z_score = - qnorm(pv/2)
    for(moment_order in 2:1){
      # Calculate the final statistics based on the threshold
      s = 0
      for (i in 1:(p-1)){
        for (j in (i+1):p){
          if ( graph[i,j]==0 & pv[i, j] < threshold){
            s = s + abs(z_score[i, j])^moment_order
          }
        }
      }
      all_rc_stat=c(all_rc_stat,s)
    }
  }
  names(all_rc_stat)=paste0(ifelse(residual_type=='PR',"P","E"), 
                            c("RC_SS","RC_SA", "RC_Z_SS", "RC_Z_SA"))
  
  return(all_rc_stat)
}

## different statistics
PRC_SS=function(x, graph,...){
  return(Stat_RC(x,graph,'PR',"exact",2,...))
}

PRC_SA=function(x, graph,...){
  return(Stat_RC(x,graph,'PR',"exact",1,...))
}

PRC_Z_SS=function(x, graph,...){
  return(Stat_RC(x,graph,'PR',"stablizing",2,...))
}
PRC_Z_SA=function(x, graph,...){
  return(Stat_RC(x,graph,'PR',"stablizing",1,...))
}


ERC_SS=function(x, graph,...){
  return(Stat_RC(x,graph,'ER',"exact",2,...))
}

ERC_SA=function(x, graph,...){
  return(Stat_RC(x,graph,'ER',"exact",1,...))
}

ERC_Z_SS=function(x, graph,...){
  return(Stat_RC(x,graph,'ER',"stablizing",2,...))
}
ERC_Z_SA=function(x, graph,...){
  return(Stat_RC(x,graph,'ER',"stablizing",1,...))
}


#Bonforrnni method by Drton, M., & Perlman, M. D. (2007). Multiple Testing and Error Control in Gaussian Graphical Model Selection as benchmark
Bonf_GoF=function(x, graph){
  p=ncol(x)
  pv=PRC_pvalues(x, graph)$pv.exact
  N0=sum(0==graph[upper.tri(graph,diag=F)])
  return(min(min(pv)*N0,1))
}


## local regression-based methods-----
local_regression_Fvalue <- function(x, graph){
  p = ncol(x)
  n = nrow(x)
  Fvalue = matrix(0, nrow = p, ncol = p)  # Initialize
  pvalue = matrix(1, nrow = p, ncol = p)
  
  
  
  if( max(colSums(graph)) + 3 > n){
    
    ## print('Null graph too large for local regression!')
  }else{
    
    # Find neighbors for each variable based on the adjacency graph
    neig_list = lapply(1:p, function(i) setdiff(which(graph[i, ] != 0), i))
    
    
    x=scale(x,center=T,scale=F)   # de-mean to avoid cbind when getting residual
    
    # Loop over all pairs of variables
    for (i in 1:p){
      neig_i=neig_list[[i]]
      di=length(neig_i)
      
      
      residual_null=get_residuals(x,i,neig_i)
      
      set_M1=setdiff(1:p,c(i,neig_i))
      
      for (j in set_M1){
        
        
        residual_larger = get_residuals(x,i,c(j,neig_i))
        numerator = sum((residual_larger-residual_null)^2)
        denominator = sum(residual_larger^2)/(n-1-1-di)
        
        
        
        # Store the F-test result
        Fvalue[i, j] = numerator/denominator
        pvalue[i, j] = pf(Fvalue[i, j],df1 = 1,df2 = n-1-1-di,lower.tail = F)
      }
    }
    return(list(Fvalue=Fvalue,pvalue=pvalue))   
  }
  
  local_regression_Fvalue_lmfunction <- function(x, graph){
    p = ncol(x)
    n = nrow(x)
    Fvalue = matrix(-1, nrow = p, ncol = p)  # Initialize
    pvalue = matrix(2, nrow = p, ncol = p)
    
    # Find neighbors for each variable based on the adjacency graph
    neig_list = lapply(1:p, function(i) setdiff(which(graph[i, ] != 0), i))
    
    # Loop over all pairs of variables
    for (i in 1:p){
      neig_i=neig_list[[i]]
      
      if(length(neig_i)>0){
        model_NE <- lm(x[, i] ~ x[, neig_i])
      }else{
        model_NE <- lm(x[, i] ~ 1)
      }
      
      set_M1=setdiff(1:p,c(i,neig_i))
      
      for (j in set_M1){
        
        # Fit the regression model for X_i onto the larger model
        model_larger <- lm(x[, i] ~ x[, c(neig_i, j)])
        
        # Perform the F-test to compare the models
        f_test_result <- anova(model_NE, model_larger)
        
        # Store the F-test result
        Fvalue[i, j] = f_test_result$F[2]
        pvalue[i, j] = f_test_result$`Pr(>F)`[2]
      }
    }
  }
  return(list(Fvalue=Fvalue,pvalue=pvalue))   
}


F_max=function(x, graph){
  Fp_stat=local_regression_Fvalue(x,graph)
  Fv=Fp_stat$Fvalue
  return(max(Fv))
}


F_sum=function(x, graph){
  Fp_stat=local_regression_Fvalue(x,graph)
  Fv=Fp_stat$Fvalue
  
  return(sum(Fv[Fv>0]))
}


F_logsum=function(x, graph){
  Fp_stat=local_regression_Fvalue(x,graph)
  Fv=Fp_stat$Fvalue
  
  return(sum(log(Fv[Fv>0])))
}

F_Z=function(x, graph,moment_order=1,threshold=0.1){  
  p=ncol(graph)
  
  Fp_stat=local_regression_Fvalue(x,graph)
  Fv=Fp_stat$Fvalue
  pv=Fp_stat$pvalue
  z_score = - qnorm(pv/2)
  
  # Calculate the final statistics based on the threshold
  FZ = 0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if ( graph[i,j]==0 & pv[i, j] < threshold){
        FZ = FZ + abs(z_score[i, j])^moment_order
      }
    }
  }
  return(FZ)
}


All_Stat_F=function(x, graph,threshold=0.1){  
  p=ncol(graph)
  
  Fp_stat=local_regression_Fvalue(x,graph)
  Fv=Fp_stat$Fvalue
  pv=Fp_stat$pvalue
  z_score = - qnorm(pv/2)
  
  # Calculate the final statistics based on the threshold
  FZ = c(0,0)
  for(moment_order in 1:2){
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        if ( graph[i,j]==0 & pv[i, j] < threshold){
          FZ[moment_order] = FZ[moment_order] + abs(z_score[i, j])^moment_order
        }
      }
    }
  }
  
  all_F_stat=c(max(Fv),sum(Fv),sum(log(Fv[Fv>0])),FZ)
  names(all_F_stat)=paste0("F_", 
                           c("max","sum", "logsum", "Z_SA","Z_SS"))
  
  return(all_F_stat)
}

# N. Verzelen and F. Villers 2009
VV_GoF=function(x,graph){
  
  Fp_stat=local_regression_Fvalue(x,graph)
  pvalues=Fp_stat$pvalue
  
  p=ncol(x)
  bonf_pvalues=c()
  for(i in 1:p){
    m_i=sum(graph[i,-i]==0)
    bonf_pvalues[i]=m_i*min(pvalues[i,])
  }
  return(min(bonf_pvalues*p,1))
}

## Likelihood ratio test ------------------
## These do not perform well

# Helper function: log likelihood of multivariate Gaussian
log_multivariate_gaussian <- function(x, covariance_matrix) {
  return(sum(dmvnorm(x, colMeans(x), covariance_matrix, log = TRUE)))
}

# Log likelihood for null model
log_likelihood_null <- function(x, graph) {
  p <- ncol(x)
  n <- nrow(x)
  
  # Initialize rho
  rho <- matrix(0, nrow = p, ncol = p)
  rho[graph == 0] <- 1e8 * n * p
  diag(rho) <- 0
  
  # Glasso
  GLS <- glasso(var(x), rho, penalize.diagonal = FALSE)
  estimated_covariance <- GLS$w
  
  return(log_multivariate_gaussian(x, estimated_covariance))
}

# Likelihood ratio using ridge regularization for full model
glr_ridge <- function(x, graph, N_CV=10) {
  p <- ncol(x)
  n <- nrow(x)
  
  # Standard deviations
  std_devs <- sqrt(diag(var(x)))
  
  # Ridge for full model
  correlation_matrix <- cor(x)
  if(N_CV>1){
    opt_lambda <- optPenalty.aLOOCV(scale(x), lambdaMin = .001, lambdaMax = 30, step = N_CV, verbose = FALSE,
                                    graph=F   )$optLambda
  }else{
    opt_lambda=p/n
  }
  inv_sigma_cor <- ridgeP(correlation_matrix, opt_lambda, type = "Alt", target = correlation_matrix)
  
  # Compute Sigma
  sigma_matrix <- sweep(sweep(solve(inv_sigma_cor), 1, std_devs, "*"), 2, std_devs, "*")
  log_likelihood_full <- log_multivariate_gaussian(x, sigma_matrix)
  
  # # Null model likelihood
  # log_likelihood_null_model <- log_likelihood_null(x, graph)
  # return(2 * (log_likelihood_full - log_likelihood_null_model))
  
  return(2 * log_likelihood_full ) 
  
}

# Likelihood ratio using glasso for full model
glr_glasso <- function(x, graph, lambda=NULL) {
  p <- ncol(x)
  n <- nrow(x)
  
  
  if(is.null(lambda)){
    lambda=sqrt(log(p) / n)
  }
  
  # Standard deviations
  std_devs <- sqrt(diag(var(x)))
  
  # Glasso for full model
  rho <- matrix(0, nrow = p, ncol = p)
  rho[graph == 0] <- lambda
  diag(rho) <- 0
  
  GLS <- glasso(cor(x), rho, penalize.diagonal = FALSE)
  sigma_matrix <- sweep(sweep(GLS$w, 1, std_devs, "*"), 2, std_devs, "*")
  
  log_likelihood_full <- log_multivariate_gaussian(x, sigma_matrix)
  
  # # Null model likelihood
  # log_likelihood_null_model <- log_likelihood_null(x, graph)
  # return(2 * (log_likelihood_full - log_likelihood_null_model))
  
  return( 2 * log_likelihood_full)
  
}

## temporary testing

if(1==2) # a piece of code for testing
{
  rm(list=ls())
  library(Rfast)
  
  n=30
  p=10
  graph= (abs(outer(1:p,1:p,FUN ='-'))<3)*1
  x=matrix(rnorm(n*p),nc=p)
  css=exchangeable_sampling(x,graph)
  
  length(css$X_copies)
  
  PRC_SS(x,graph,threshold=0.5)
  ERC_SS(x,graph,threshold=0.5)
  dp_bonf(x,graph)
  
  Stat_LR(x,graph,threshold = 0.5)
  glr_glasso(x,graph)
  glr_ridge(x,graph)
}




