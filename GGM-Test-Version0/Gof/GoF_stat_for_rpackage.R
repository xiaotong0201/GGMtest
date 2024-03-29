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
          
        }else {
          pv.exact[i, j] = 1
        }
        # Store symmetric values
        pv.exact[j, i] = pv.exact[i, j]
      }
    }
  }
  return(list(pv.exact=pv.exact))
}

PRC_SS=function(x, graph,moment_order=2,threshold=0.1){
  p=ncol(x)
  
  pv_list=PRC_pvalues(x, graph)
  pv=pv_list$pv.exact
  
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



# Compute the Efficient Residual Correlation (ERC) p-values
ERC_pvalues <- function(x, graph){
  p = ncol(x)
  n = nrow(x)
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
          pv.stab[i, j] = pvalue_stablizing_cor(res_cor,df)
          
          
        } else {
          pv.stab[i,j] = 1
        }
        # Store symmetric values
        pv.stab[j, i] = pv.stab[i, j] 
      }
    }
  }
  return(list(pv.stab=pv.stab))   
}

ERC_Z_SS=function(x, graph,moment_order=2,
                  threshold=0.1){
  p=ncol(x)
  pv_list=ERC_pvalues(x,graph)
  pv=pv_list$pv.stab
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


# F_sum

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

F_sum=function(x, graph){
  Fp_stat=local_regression_Fvalue(x,graph)
  Fv=Fp_stat$Fvalue
  return(sum(Fv[Fv>0]))
}


# Likelihood ratio using glasso for full model

# Helper function: log likelihood of multivariate Gaussian
log_multivariate_gaussian <- function(x, covariance_matrix) {
  return(sum(dmvnorm(x, colMeans(x), covariance_matrix, log = TRUE)))
}

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

  return( 2 * log_likelihood_full)
  
}
