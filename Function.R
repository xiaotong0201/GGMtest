

# Function for Algorithm 1 and 2 --------------------------------------------------------

################ Sampling Exchangeable Copies ################

### Residual Rotation
### Algorithm 1

# Generates a uniform vector orthogonal to the given matrix A.

ortho.vec = function(A) {
  # Perform QR decomposition to obtain orthonormal columns
  qr.dec = qr(cbind(A, rnorm(nrow(A))))
  # Generate the orthogonal vector
  return(qr.Q(qr.dec)[, 1+ncol(A)] * sample(c(-1, 1), 1))
}

# Rotates one coordinate.

rotate = function(x, i, graph) {
  # Find neighbors in the graph for the i-th coordinate
  neighbors = setdiff(which(graph[i, ] != 0), i)
  n = nrow(x)
  vec.1 = matrix(1, nr=n, nc=1)
  new.xi = x[, i]
  sample.mean = mean(x[, i])
  
  
  l=length(neighbors)
  if(l == 0) {
    # If no neighbors, rotate based on the sample mean
    new.xi = sample.mean + ortho.vec(vec.1) * sqrt(sum((x[, i] - sample.mean)^2))
  } else if(n >=  2+l){

    lr.temp = lm(x[, i]~ x[, neighbors] )
    new.xi = lr.temp$fitted.values +  sqrt(sum(lr.temp$residuals^2)) * ortho.vec(cbind(vec.1, x[, neighbors])) 
   
  }
  
  return(new.xi)
}

# Runs one round of the algorithm

run = function(x, orders, graph) {
  newx = x
  for(i in orders) {
    newx[, i] = rotate(newx, i, graph)
  }
  return(newx)
}

### Sampling Exchangeable Copies from Markov Chains
### Algorithm 2

# Sampling exchangeable copies for a graph

exchangeable_sampling <- function(X, graph, M=100, L=1, 
                                  I=seq(1:ncol(X)), bReduced = F) {
  
  # Step 1: Start from X and run the previously defined 'run' function 
  # according to the order of I for L times to generate X_hub.
  X_hub <- X
  for (l in 1:L) {
    X_hub <- run(X_hub, I, graph)
  }
  
  # Initialize an empty list to store the M copies
  X_copies <- vector("list", M)
  
  # Step 2: For each copy, start from X_hub and run the 'run' function
  # according to the reversed order of I for L times to generate each X_tilde.
  for (m in 1:M) {
    X_tilde <- X_hub
    for (l in 1:L) {
      X_tilde <- run(X_tilde, rev(I), graph)
    }
    if(bReduced){ 
      X_copies[[m]] <- X_tilde[,I]  # no need to store others for CRT
    }else{
      X_copies[[m]] <- X_tilde
    }
  }
  
  # Output: X_copies containing the M transformed matrices
  return(list(X_hub=X_hub, X_copies=X_copies))
}

# Compute the p-value based on observed and sampled test statistics
ComputePValue <- function(T0, Ttilde, type =  "One-sided", #c('One-sided', 'Two-sided')
                          randomize = TRUE) {
  
  
  if(is.matrix(Ttilde)){
    cntGreater <- rowSums(Ttilde > T0)
    cntSmaller <- rowSums(Ttilde < T0)
    cntEqual   <- rowSums(Ttilde == T0) + 1
    M <-  ncol(Ttilde)
  }else{
    
    cntGreater <- sum(Ttilde > T0)
    cntSmaller <- sum(Ttilde < T0)
    cntEqual   <- sum(Ttilde == T0) + 1
    M <- length(Ttilde)
  }
  
  # Random sampling for tie-breaking
  S <- ifelse(randomize, 
              sapply(cntEqual,function(x)sample(x,1)) # cntEqual may be a vector
              , cntEqual)
  
  # Calculate p-value
  if (type == "One-sided") {
    pvalue <- (cntGreater + S) / (1 + M)
  } else {
    pvalue <- 2 * ( pmin(cntGreater, cntSmaller) + S) / (1 + M)
  }
  
  return(pvalue)
}



# Function for GoF --------------------------------------------------------

################ Goodness-of-fit Test for GGMs ################

### GoF Statistics

## Residual correlation based methods 

# Helper function to get residuals
# Assume x has been de-mean
get_residuals <- function(x,i, neighbors) {
  l=length(neighbors) 
  n=nrow(x)
  if (l== 0) {
    return(x[,i])
  }else if(1+l>= n){
    return( x[,i]*0 )
  }else {
    model = lmfit(x[, neighbors], x[,i])
    return(model$residuals)
  }
}

# Assume intercept is needed
get_residuals_intercept <- function(x,i, neighbors) {
  if (length(neighbors) == 0) {
    return(x[,i] - mean(x[,i]))
  }else if(length(neighbors)+1>n){
    return( x[,i]*0 )
  } else {
    model = lmfit(cbind(1,x[, neighbors]), x[,i])
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
  pv.stab = matrix(1, nrow = p, ncol = p)
  
  
  
  x=scale(x,center=T,scale=F)   # de-mean to avoid cbind when getting residual
  
  # Find neighbors for each variable based on the adjacency graph
  neig_list = lapply(1:p, function(i) setdiff(which(graph[i, ] != 0), i))
  
  # Loop over all pairs of variables
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if (graph[i, j] == 0){  # if i and j are not adjacent
        
        # Find the union of neighbors for i and j
        neig = union(neig_list[[i]], neig_list[[j]])
        df = n - length(neig) - 2
        
        if(df>= 1){
          # Calculate residuals
          res_i = get_residuals(x, i, neig)
          res_j = get_residuals(x, j, neig)
          
          # Compute residual correlation and test statistics
          res_cor = cor(res_i,res_j)
          
          
          # Compute p-value 
          t_stat = sqrt(df) * res_cor / sqrt(1 - res_cor^2)
          pv.exact[i, j] = 2 * pt(abs(t_stat), df, lower.tail=FALSE)
          pv.stab[i,j]=pvalue_stablizing_cor(res_cor,df)
        }
        
        
      }
      # Store symmetric values
      pv.exact[j, i] = pv.exact[i, j]
      pv.stab[j, i] = pv.stab[i, j] 
      
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
        
        
        df=n-2-min(size_neig[i],size_neig[j])
        if(df>=1){
          
          res_cor = cor(res_i,res_j)
          pv.exact[i, j] = pbeta(abs(res_cor),shape1=1,shape2=df,lower.tail = F)
          pv.stab[i, j] = pvalue_stablizing_cor(res_cor,df)
        }
        
      } 
      # Store symmetric values
      pv.exact[j, i] = pv.exact[i, j]
      pv.stab[j, i] = pv.stab[i, j] 
    }
  }
  
  return(list(pv.exact=pv.exact,pv.stab=pv.stab))   
}

# compute one statistic

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

## Bonforrnni method by Drton, M., & Perlman, M. D. (2007). Multiple Testing and Error Control in Gaussian Graphical Model Selection as benchmark

Bonf_GoF=function(x, graph){
  p=ncol(x)
  pv=PRC_pvalues(x, graph)$pv.exact
  N0=sum(0==graph[upper.tri(graph,diag=F)])
  return(min(min(pv)*N0,1))
}

## local regression-based methods

local_regression_Fvalue <- function(x, graph){
  p = ncol(x)
  n = nrow(x)
  Fvalue = matrix(0, nrow = p, ncol = p)  # Initialize
  pvalue = matrix(1, nrow = p, ncol = p)
  
  
  
  
  # Find neighbors for each variable based on the adjacency graph
  neig_list = lapply(1:p, function(i) setdiff(which(graph[i, ] != 0), i))
  
  
  x=scale(x,center=T,scale=F)   # de-mean to avoid cbind when getting residual
  
  # Loop over all pairs of variables
  for (i in 1:p){
    neig_i=neig_list[[i]]
    di=length(neig_i)
    
    if(n >= 3+di){
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
  num_fullnode=0  # Ignore nodes that have full connections
  bonf_pvalues=rep(1,p)
  for(i in 1:p){
    m_i=sum(graph[i,-i]==0)
    if(m_i>0){
      bonf_pvalues[i]=m_i*min(pvalues[i,])
    }else{
      num_fullnode=num_fullnode+1
    }
  }
  return(min(bonf_pvalues*(p-num_fullnode),1))
}

## Likelihood ratio test

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

### Conducts a goodness-of-fit test based on exchangeable sampling
### Algorithm 3

CSSGoF <- function(X, graph, X_copies = NULL, 
                   Fun_Stat, type='One-sided', randomize=T,M=100,L=1) {
  
  # Generate exchangeable samples if not provided
  if (is.null(X_copies)) {
    CSS_samples <- exchangeable_sampling(X, graph, M=M, L=L)
    X_copies <- CSS_samples$X_copies
  }
  
  T0 <- Fun_Stat(X, graph)
  Ttilde <- simplify2array(lapply(X_copies, Fun_Stat, graph))
  
  pvalue=ComputePValue(T0, Ttilde, type=type, randomize = randomize)
  return(list(pvalue=pvalue, T0=T0, Ttilde=Ttilde))
}

# Function to simulate the data
simulate_data <- function(Sigma, n) {
  X <- mvrnorm(n, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
  return(X)
}

GGM =function(popul_param,
              CSS_param=list(M=100,L=3),
              epo, 
              List_CSS_stat=NULL,
              CSS=NULL,dir_prefix=''
){
  
  
  if(is.null(List_CSS_stat)==T){
    
    List_CSS_stat=c(
      "PRC_SS", "PRC_SA", "PRC_Z_SS", "PRC_Z_SA",
      "ERC_SS", "ERC_SA", "ERC_Z_SS", "ERC_Z_SA",
      "F_max", "F_sum","F_logsum", "F_Z_SA","F_Z_SS", 
      "glr_ridge", "glr_glasso"
    )
    
  }
  
  p <- popul_param$P
  n <- popul_param$N
  model <- popul_param$model
  modelparam <- popul_param$modelParam
  
  #print(paste(n,p,model,paste0(unlist(modelparam),collapse = ','),sep=","))
  
  path_prefix=paste0(dir_prefix,'DataStorage-',format(Sys.time(), "%b-%e-%Y"),'/')
  #     # Create the folder if the folder not exist
  if (!file.exists(path_prefix)) {
    dir.create(path_prefix)
  }
  file_prefix=paste0(c(rbind(names(popul_param)[1:3],unlist(popul_param[1:3]))),collapse ="_")
  
  file_prefix=paste0(path_prefix,file_prefix,'_',
                     paste0(c(rbind(names(modelparam),unlist(modelparam))),collapse ="_")
  )
  
  
  set.seed(epo) 
  
  
  # Setting up the population and the null graph
  if (model=="Band_Band"){
    true_graph=graph_band(p,modelparam$K,modelparam$s)
    null_graph=graph_band(p,modelparam$K0,modelparam$s) 
  }else if(model == "Hub"){
    true_graph=graph_hub(p,modelparam$hubsize,0,modelparam$noise)
    null_graph=graph_hub(p,modelparam$hubsize,modelparam$Q,modelparam$noise)
  }else if(model == "ER"){
    true_graph=graph_ER(p,modelparam$s,modelparam$Q,epo)
    null_graph=graph_ER_null(p,modelparam$Q,modelparam$Q0,true_graph)
  }else if(model== "Verzelen"){
    true_graph=graph_Verzelen(p,modelparam$E)
    null_graph=graph_Verzelen_null(p,modelparam$Q,modelparam$E,true_graph$inv.sigma)
  }else{
    cat('Error: Setting up the model.')
    return(NULL)
  }
  
  sigma=true_graph$sigma
  inv.sigma=true_graph$inv.sigma
  G0=prec_to_adj(null_graph$inv.sigma)
  
  if(min(eigen(inv.sigma)$values)<0){
    print('Model Error!')
    return(NULL)
  }
  
  # Simulate the data
  x <- simulate_data(sigma, n)
  saveRDS(x, file = paste0(file_prefix,"_n=",n,"_p=",p,"_x_epo=",epo,".rds") )
  
  ### Inference begins
  
  ## Verzelen et al 
  VV_pvalue <- VV_GoF(x,G0)
  ## Bonferroni
  DP_pvalue <- Bonf_GoF(x,G0)
  
  
  ### CSS sampling
  if(is.null(CSS)){
    start.time <- Sys.time()
    #cat('CSS begin;')
    CSS=exchangeable_sampling(x,G0, M=CSS_param$M, L=CSS_param$L)
    time.taken <- Sys.time() - start.time
    #cat('Time used: ', time.taken,'s;')
    #cat('CSS end;\n')
    saveRDS(CSS, file = paste0(file_prefix,"_n=",n,"_p=",p,"_CSS_epo=",epo,".rds") )
  }
  
  #calculate statistics
  #cat('CSS testing begin;')
  CSS_test_result=list()
  
  
  i=1
  while(i<= length(List_CSS_stat)){
    f=List_CSS_stat[i]
    start.time <- Sys.time()
    if(substr(f,2,3)=="RC"){
      temp_fn=function(x,g)All_Stat_RC(x,g,residual_type=substr(f,1,2))
      bundle_results=CSSGoF(x,G0, CSS$X_copies, Fun_Stat = temp_fn)
      for(j in 1:length(bundle_results$pvalue)){
        CSS_test_result[[i+j-1]] = list(
          pvalue=bundle_results$pvalue[j],
          T0=bundle_results$T0[j],
          Ttilde=bundle_results$Ttilde[j,]
        )
      }
      i=i+3
    }else if(startsWith(f,"F_")){
      bundle_results=CSSGoF(x,G0, CSS$X_copies, Fun_Stat = All_Stat_F)
      for(j in 1:length(bundle_results$pvalue)){
        CSS_test_result[[i+j-1]] = list(
          pvalue=bundle_results$pvalue[j],
          T0=bundle_results$T0[j],
          Ttilde=bundle_results$Ttilde[j,]
        )
      }
      i=i+4
      
    }else{
      CSS_test_result[[i]] = CSSGoF(x,G0, CSS$X_copies, Fun_Stat = get(f))
    }
    
    
    time.taken <- Sys.time() - start.time
    #cat('Time used: ', time.taken,'s;')
    CSS_test_result[[i]]$time.taken=time.taken
    
    i=i+1
  }
  
  #cat('CSS testing end;\n')
  saveRDS(CSS_test_result, file = paste0(file_prefix,"_n=",n,"_p=",p,"_TestStat_epo=",epo,".rds") )
  
  CSS_test_pvalues = unlist(lapply(CSS_test_result, function(a)a$pvalue))
  names(CSS_test_pvalues)=List_CSS_stat
  
  # Wrap results in a list
  final_result <- list(
    'VV_pvalue' = VV_pvalue,
    'DP_pvalue' = DP_pvalue,
    'CSS_test_pvalues' = CSS_test_pvalues
  )
  return(final_result)
}     

# Experiment

run_experiments_epo <- function(param_template, epo, M,L,List_CSS_stat=NULL,
                                dir_prefix='') {
  final_results <- list()
  lFR=0
  for (N in param_template$N) {
    for (P in param_template$P) {
      # Create a list to hold the possible combinations of model parameters
      combos <- expand.grid(param_template$modelParam)
      
      # Loop over each combination of model parameters
      for (i in 1:nrow(combos)) {
        
        modelParam = as.list(combos[i, ])
        names(modelParam)=names(param_template$modelParam)
        
        param_part <- list(model=param_template$model, 
                           N=N,P = P, modelParam=modelParam)
        # Call GGM function and store the result
        final_results[[lFR+1]] = GGM(param_part, 
                                     CSS_param=list(M=M,L=L),
                                     epo=epo, List_CSS_stat=List_CSS_stat,
                                     dir_prefix=dir_prefix)
        lFR=lFR+1
      }
    }
  }
  return(final_results)
}

# read data
read_experiments_epo <- function(param_template, output_list,sig.lv=0.05) {
  MaxEpo=length(output_list)
  
  summary_result <- list()
  lFR=0
  
  for (ind.N in 1:length(param_template$N)) {
    N=param_template$N[ind.N]
    for (ind.P in 1:length(param_template$P)) {
      P=param_template$P[ind.P]
      
      combos <- expand.grid(param_template$modelParam)
      
      # Loop over each combination of model parameters
      for (ind.setup in 1:nrow(combos)) {
        
        
        np.id=(ind.N-1)*length(param_template$P)+ind.P
        exp.id=ind.setup+nrow(combos)*(np.id-1)
        
        pvalue.table=c()
        for(epo in 1:MaxEpo){
          pvalue.table= rbind(pvalue.table, c(unlist(output_list[[epo]][[exp.id]]) ))
        }
        
        colnames(pvalue.table)=c("VV","DP",names(output_list[[1]][[1]]$CSS_test_pvalues))
        
        # colMeans(pvalue.table)
        power=colMeans(pvalue.table < sig.lv)
        SE=apply(pvalue.table < sig.lv,2,function(v)sqrt(var(v)/length(v)))
        
        # print(power)
        
        modelParam = as.list(combos[ind.setup, ])
        names(modelParam)=names(param_template$modelParam)
        
        summary_result[[lFR+1]]  <- 
          list(model=param_template$model, 
               N=N,P = P, modelParam=modelParam,
               pvalue.table=pvalue.table,
               power=power,SE=SE
          )
        lFR=lFR+1
        
      }
    }
  }
  return(summary_result)
}


extract_metric=function(metric_name, sig.lv=0.05){
  if(metric_name=='ospt'){ # one sample proportional test
    one.sample.prop.test=sapply(summ, 
                                function(res){
                                  apply(res$pvalue.table,2, 
                                        function(vv){
                                          
                                          prop.test( sum(vv < sig.lv), length(vv), p=0.05, alternative='greater',correct=F)$p.value
                                          
                                        }
                                  )
                                })
    new_table=sapply(summ,function(x)(c(x$N,x$P,unlist(x$modelParam),apply(x$pvalue.table,2, 
                                                                           function(vv){
                                                                             prop.test( sum(vv < sig.lv), length(vv), p=0.05,alternative='greater',correct=F)$p.value
                                                                           })
    )))
  }else{
    new_table=sapply(summ,function(x)(c(x$N,x$P,unlist(x$modelParam),x[[metric_name]])))
  }
  rownames(new_table)[1:2]=c('n','p')
  new_table=t(new_table)
  new_table=as.data.frame(new_table)
  return(new_table)
}

################ Graph Setting ################

# Helper function 
#  convert a precision matrix to an adjacency matrix for a graph
prec_to_adj <- function(prec_mat,  threshold = 1e-7) {
  adj_mat= 1*(abs(prec_mat) > threshold)
  diag(adj_mat)=0
  return(adj_mat)
}

# Band model
# K times s should be smaller than 2
graph_band=function(p,k,s){
  if(k*abs(s)>=2){
    print('Model Error (Band)')
    return(NULL)
  }
  inv.sigma <- diag(p)
  for (i in 1:p){
    for (j in 1:p){
      if(abs(i - j) >= 1 && abs(i - j) <= k){
        inv.sigma[i,j]=s
        inv.sigma[j,i]=s
      }
    }
  }
  return(list(sigma = solve(inv.sigma), inv.sigma = inv.sigma))
}

# Hub Model
# delete an edge with probability q
graph_hub=function(p,hubsize,q=0,noise=0.5){
  
  Omega1 <- matrix(0, nrow = p, ncol = p)
  for (k in 1:(p/hubsize)) {
    i <- hubsize * (k - 1) + 1
    for( j in (i+1):(i+hubsize-1) ){
      if(rbinom(1,1,q)==0){  
        Omega1[i,j]=1
        Omega1[j,i]=1
      }
    }
  }
  diag(Omega1)=colSums(abs(Omega1)) + noise
  return(list(sigma = solve(Omega1), inv.sigma = Omega1))
}

# Erdos–Renyi model
# Edge connected with probability q
graph_ER=function(p,s,q,epo){
  set.seed(epo)
  inv.sigma <- matrix(0,p,p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      signal=runif(1, min = s/2, max = 3*s/2)
      if (runif(1) < q){
        inv.sigma[i,j]= signal
        inv.sigma[j,i]= signal
      }
    }
  }
  # Make it positive definite as per LaTeX description
  lambda <- min(eigen(inv.sigma)$values)
  inv.sigma = inv.sigma + (abs(lambda) + 0.05) * diag(p)
  return(list(sigma = solve(inv.sigma), inv.sigma = inv.sigma))
}

# marginally, the edge connect with probability q0
graph_ER_null=function(p,q,q0,er_graph){
  inv.sigma <- er_graph$inv.sigma
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if (runif(1) > q0/q){  # delete w.p. 1- q0/q
        inv.sigma[j,i]=inv.sigma[i,j] = 0
      }
    }
  }
  return(list(inv.sigma=inv.sigma))
}

#Reproduce the example in N. Verzelen and F. Villers 2009

graph_Verzelen= function(p, eta) {
  
  
  num_edges <- floor(eta * p * (p - 1) / 2) # tricky! 2023-11-01
  adjacency_matrix <- matrix(0, nrow = p, ncol = p)
  
  edges <- sample(which(lower.tri(adjacency_matrix,diag = F)), num_edges)
  
  adjacency_matrix[edges] <- runif(num_edges, -1, 1)
  
  adjacency_matrix=adjacency_matrix+t(adjacency_matrix)
  
  Pi = generate_partial_corr_matrix(adjacency_matrix) 
  inv.sigma = -Pi; diag(inv.sigma)=1 
  
  ## standardization
  Sigma=solve(inv.sigma)
  diag_s=diag(Sigma)
  Sigma=Sigma/ sqrt(diag_s %o% diag_s)
  ##
  
  
  return(list(sigma = Sigma, inv.sigma = inv.sigma) )
}

generate_partial_corr_matrix= function(U) {
  diag_values <- apply(abs(U), 2, sum) + 1e-4
  diag(U) <- diag_values
  Pi <- U / sqrt(diag_values %o% diag_values)
  return(Pi)
}

graph_Verzelen_null=function(p,q,eta,prec){
  g_s= 1*(prec!=0)
  
  num_edges_to_delete <- ceiling(q * p * (p - 1) / 2 * eta)  # tricky! 2023-11-01
  
  gs_ind=which( (g_s != 0)&lower.tri(g_s,diag=F), arr.ind = TRUE)
  edges_to_delete <- sample(nrow(gs_ind), min(nrow(gs_ind), num_edges_to_delete))
  g_minus_q_s <- g_s
  g_minus_q_s[gs_ind[edges_to_delete,]] <- 0
  g_minus_q_s = g_minus_q_s* (t(g_minus_q_s!=0))
  return(list(inv.sigma=g_minus_q_s))
  
}

# Function for CRT --------------------------------------------------------

################ Conditional Randomization Test ################

# Helper function
# Dependencies: glmnet 
REF_DS_inf <- function(x, y, family, lasso_est) {
  nn <- length(y)
  X <- cbind(rep(1, nrow(x)), x)
  if(ncol(X) != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }
  if(family == "binomial") {
    mu <- as.vector(exp(X%*%lasso_est)/(1+exp(X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu*(1-mu))%*%X/nn
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
  } else {
    stop("Input family is not supported.")
  }
  
  theta_inv <- solve(neg_ddloglik_glmnet)
  b_hat_inv <- as.vector(lasso_est - theta_inv%*%neg_dloglik_glmnet)
  se_inv <- sqrt(diag(theta_inv))/sqrt(nn)
  pval_inv <- 2*pnorm(abs(b_hat_inv/se_inv), lower.tail=F)
  
  return(list(est=b_hat_inv, se=se_inv, pvalue=pval_inv, theta=theta_inv))
}


ORIG_DS_inf <- function(x, y, family, lasso_est,beta_glmnet, nfold=5, n_lambda=100, lambda_ratio=0.005) {
  nn <- length(y)
  pp <- ncol(x)
  X <- cbind(rep(1, nrow(x)), x)
  if(ncol(X) != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }
  if(family == "binomial") {
    mu <- as.vector(exp(X%*%lasso_est)/(1+exp(X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu*(1-mu))%*%X/nn
    C_glmnet <- sqrt(diag(mu*(1-mu))/nn)%*%X
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
    C_glmnet <- sqrt(diag(mu)/nn)%*%X
  } else {
    stop("Input family is not supported.")
  }
  
  theta_glmnet <- diag(pp+1)
  tau_glmnet <- rep(NA, pp+1)
  for(j in 1:(pp+1)) { # for: nodewise lasso 
    current_x <- sqrt(nn)*C_glmnet[,-j]
    current_y <- sqrt(nn)*as.vector(C_glmnet[,j])
    lam_max <- max(abs(t(current_x)%*%current_y)/nn)
    lam_min <- lam_max*lambda_ratio
    lam_seq <- exp(seq(from=log(lam_max), to=log(lam_min), length.out=n_lambda))
    gamma_j_glmnet <- cv.glmnet(x=current_x, y=current_y,
                                family="gaussian", alpha=1, standardize=F, intercept=F,
                                nfolds=nfold, lambda=lam_seq)
    gamma_j_glmnet <- as.vector(glmnet(x=sqrt(n)*C_glmnet[,-j], y=sqrt(n)*as.vector(C_glmnet[,j]),
                                       family="gaussian", alpha=1, standardize=F, intercept=F,
                                       lambda=gamma_j_glmnet$lambda.min)$beta)
    theta_glmnet[j,-j] <- (-1)*t(gamma_j_glmnet)
    tau_glmnet[j] <- as.numeric(neg_ddloglik_glmnet[j,j]-neg_ddloglik_glmnet[j,-j]%*%gamma_j_glmnet)
  } # end for: nodewise lasso
  theta_glmnet <- diag(1/tau_glmnet)%*%theta_glmnet 
  
  b_hat_nw <- as.vector(beta_glmnet - theta_glmnet%*%neg_dloglik_glmnet)
  se_nw <- sqrt(diag(theta_glmnet%*%neg_ddloglik_glmnet%*%t(theta_glmnet)))/sqrt(nn)
  pval_nw <- 2*pnorm(abs(b_hat_nw/se_nw), lower.tail=F)
  
  return(list(est=b_hat_nw, se=se_nw, pvalue=pval_nw, theta=theta_glmnet))
}

### CRT-stat

## For direct computations 

SST_GLM=function(X,Y,setT=NULL,model_type="gaussian",...)
{
  if(is.null(setT)){setT=1:ncol(X)}
  glm_fit <- glm( Y ~ X, family = model_type)
  T_stat=summary(glm_fit)$coefficients[1+setT,3]
  return(sum(T_stat^2))
}


Dev_GLM=function(X,Y, setT=NULL,model_type="gaussian",...)
{
  glm_fit <- glm( Y~X, family = model_type)
  return( - glm_fit$deviance)
}


F_GLM=function(X,Y, setT=NULL,model_type="gaussian",...)
{
  glm_fit <- glm( Y~X, family = model_type)
  glm_fit_null <- glm( Y~X[,-setT,drop=F], family = model_type)
  
  return(1/anova(glm_fit_null,glm_fit,test = 'F')[2,6] )
}

library(randomForest)

RF=function(X,Y, setT=NULL,model_type='regression', RF.num.trees=100,...)
{
  
  if(model_type=="classification"){
    Y=as.factor(Y)
  }
  
  if(is.null(setT)){setT=1:ncol(X)}
  
  rf <- randomForest(X, Y, ntree=RF.num.trees, importance=TRUE) 
  if(model_type=="regression"){
    feature_importance=rf$importance[,1]
  }else{
    feature_importance=rf$importance[, ncol(rf$importance)-1]
  }
  return(sum(pmax(feature_importance[setT],0)))
}

## For two-step computation 

GLM_NULL=function(X,Y, null_model_type="gaussian",...)
{
  glm_fit <- glm(Y~X, family = null_model_type)
  return(list(fitted.value=glm_fit$fitted.values))
}

Lasso_CV_Null=function(X,Y,null_model_type="gaussian",nfolds=5,k=0,lambda_method="lambda.1se",...)
{
  
  n=length(Y)
  cv_lasso <- cv.glmnet(X, Y,family = null_model_type, nfolds = nfolds)
  lasso = glmnet(X, Y, family = null_model_type)
  fitted.value=predict(lasso,X,s=cv_lasso[[lambda_method]],type="response")
  if(k<0){
    imp.var=NULL
  }else{
    
    beta_hat = coef(lasso,s=cv_lasso[[lambda_method]])
    
    if (null_model_type == "multinomial") {
      beta_hat_sq = rowSums(simplify2array(sapply(beta_hat, function(v)as.numeric(v^2))))
      imp.var <- which(beta_hat_sq[-1] != 0) # return the important variable pick by CV
    } else {
      imp.var <- which(beta_hat[-1] != 0) # return the important variable pick by CV
    }
    
    if(k>0 | length(imp.var) > n/2 ){ # k is nonzero or cv pick too many
      
      find_first_nonzero <- function(beta) {
        apply(beta, 1, function(x) match(TRUE, abs(x) > 0))
      }
      if (null_model_type == "multinomial") {
        indices <- apply(sapply(lasso$beta, find_first_nonzero), 1, min)
      } else {
        indices <- find_first_nonzero(lasso$beta)
      }
      imp.var <- order(indices)[1:  max(k , n/2) ]
    }
  }
  return(list(fitted.value=fitted.value,imp.var=imp.var))
}



RF_Null=function(X,Y, null_model_type="regression", RF.num.trees=100, k = 0.8,... )
{
  if(null_model_type=="classification"){
    Y=as.factor(Y)
  }
  
  rf <- randomForest(X, Y, ntree=RF.num.trees, importance=TRUE) 
  fitted.value=rf$predicted
  
  
  if(k<= 0){
    imp.var =NULL 
  }
  else{
    
    if(null_model_type=="regression"){
      feature_importance=rf$importance[,1]
    }else{
      feature_importance=rf$importance[, ncol(rf$importance)-1]
    }
    sort.imp=sort(feature_importance,decreasing = T, index.return = T)
    
    if(k<1){
      cum.imp=cumsum(sort.imp$x)
      if(cum.imp[1]>0){
        k=min(match(T,cum.imp>= k* max(cum.imp)))
      }
    }
    if(k>=1){
      imp.var=sort.imp$ix[1:k]
    }else{
      imp.var=c()
    }
  }
  return(list(fitted.value=fitted.value,imp.var=imp.var))
}

## Baseline Functions


bonfTest <- function(pv){
  return(min(min(pv)* length(pv), 1))
}

ANOVA_GLM <- function(X,Y,setT,model_type){
  red = glm(Y ~ X[, -setT,drop=F],family = model_type)
  full = glm(Y ~ X, family = model_type)
  
  if(model_type=='gaussian'){
    anova_result=anova(red, full, test='F')
    pv =  anova_result[2,6]
  }else{
    anova_result=anova(red, full, test='Chisq')
    pv =  anova_result[2,5]  
  }
  
  return(pv)
}

library(VGAM)

Chisq_Multi=function(X,Y,setT,...){
  
  
  reduced = vglm(Y ~ X[, -setT,drop=F],family = multinomial)
  full = vglm(Y ~ X, family = multinomial)
  
  
  lrt=lrtest(full, reduced)
  
  return(lrt@Body[2,5])
}

library(hdi)
DL_hdi <- function(X,Y,setT,
                   model_type="gaussian" # can be binomial or gaussian
){
  pv = bonfTest(lasso.proj(X, Y, family = model_type,
                           multiplecorr.method = "holm")$pval[setT]
  )
  return(pv)
}

library(SIHR)
#Cai, Guo, Ma 2021
CGM=function(X,Y,setT, model_type="linear"){
  p=ncol(X)
  fit=LF(X,Y, diag(p)[,setT,drop=F], 
         model = model_type,
         intercept = TRUE,
         intercept.loading = FALSE
  )
  individual_pvalue=2*pnorm(abs(fit$est.debias.vec/fit$se.vec),lower=F)
  return(bonfTest(individual_pvalue))
}


source("DBL_GLMs_functions.R")
DL_Xia <- function(X,Y,setT,model_type,
                   debias_method="REF" # "ORIG"
){
  if(!model_type%in%c("binomial","poisson")){
    stop("Only support binomial and poisson. ")
  }
  if(!debias_method%in% c("REF","ORIG")){
    stop("Unknown debiasing method. ")
  }
  scale_X=scale(X)
  cv_glm = cv.glmnet(scale_X, Y, family=model_type, alpha=1, standardize=F, intercept=T)
  beta = as.vector(coef(glmnet(scale_X, Y, family=model_type, alpha=1, standardize=F, intercept=T), s=cv_glm$lambda.min))
  
  debias_lasso=switch(debias_method, 
                      REF  = REF_DS_inf(X, Y, family=model_type, lasso_est=beta),
                      ORIG=ORIG_DS_inf(X, Y, family=model_type, lasso_est=beta)
  )
  pv = bonfTest(debias_lasso$pvalue[setT+1])
  return(pv)
}

Individual_Ftest <- function(X,Y,setT){
  p=ncol(X)
  
  setS <- setdiff(1:p,setT)
  if(length(setS)>0){
    model_null <- lm(Y ~ X[,setS,drop=F])
  }else{
    model_null <- lm(Y ~ 1)
  }
  
  pvalues=c()
  for(j in setT){
    model_j <- lm(Y ~ X[, c(setS, j),drop=F])
    f_test_result <- anova(model_null, model_j)
    pvalues =c(pvalues, f_test_result$`Pr(>F)`[2])
  }
  return(bonfTest(pvalues))
}
source('dCRT_function.R')

dCRT <- function(X,Y,setT, model_type){
  dCRT_res = d_CRT(X,Y,FDR = 0.1, model = model_type,CDR='No',candidate_set=setT)
  return(bonfTest(dCRT_res$p_values))   # a bug fixed on 2023-11-07 
}


### Setting the test statistics for each model

# Linear regression

Stat_Comb_Linear_Dev=list(Fun_Stat = Dev_GLM, model_type="gaussian", 
                          Stat_Procedure = "Direct")  

Stat_Comb_Linear_SST=list(Fun_Stat = SST_GLM, model_type="gaussian", 
                          Stat_Procedure = "Direct")  

List_Stat_Comb_Linear=c(
  'Stat_Comb_Linear_SST'    ## LM-SST
  ,'Stat_Comb_Linear_Dev'   ## LM-SSR
)

Stat_Comb_Lasso_Res_Linear_SST=list(Fun_Stat = SST_GLM, model_type="gaussian",
                                    Fun_FitNull=Lasso_CV_Null,null_model_type="gaussian", 
                                    Stat_Procedure = "Residual")
Stat_Comb_Lasso_Res_Linear_Dev=list(Fun_Stat = Dev_GLM, model_type="gaussian",
                                    Fun_FitNull=Lasso_CV_Null,null_model_type="gaussian", 
                                    Stat_Procedure = "Residual")


List_Stat_Comb_Linear_HD=c(
  'Stat_Comb_Lasso_Res_Linear_SST',  ## LM-L1-R-SST
  "Stat_Comb_Lasso_Res_Linear_Dev"   ## LM-L1-R-SSR
)  

# Non-Linear regression

Stat_Comb_RF_R_IS=list(Fun_Stat = RF, model_type="regression", 
                       Stat_Procedure = "Direct")  


Stat_Comb_RF_R_Res_RF=list(Fun_Stat = RF, model_type="regression", 
                           Fun_FitNull = RF_Null, null_model_type="regression", 
                           Stat_Procedure = "Residual") 


Stat_Comb_RF_R_Distill_RF=list(Fun_Stat = RF, model_type="regression", 
                               Fun_FitNull = RF_Null, null_model_type="regression", 
                               Stat_Procedure = "Distill") 
List_Stat_Comb_Nonlinear=c(
  'Stat_Comb_RF_R_IS',           ## RF
  "Stat_Comb_RF_R_Distill_RF",   ## RF-RR    
  "Stat_Comb_RF_R_Res_RF"        ## RF-D
)  

List_Stat_Comb_Nonlinear_HD=c(
  List_Stat_Comb_Nonlinear
) 

# Logistic regression

Stat_Comb_Logistic_Dev=list(Fun_Stat = Dev_GLM, model_type="binomial", 
                            Stat_Procedure = "Direct")  


Stat_Comb_Logistic_Distill_Dev=list(Fun_Stat = Dev_GLM, model_type="binomial", 
                                    Fun_FitNull=Lasso_CV_Null, null_model_type="binomial", 
                                    Stat_Procedure = "Distill")

Stat_Comb_Sel_Logistic_SST=list(Fun_Stat = SST_GLM, model_type="gaussian",
                                Fun_FitNull=Lasso_CV_Null,null_model_type="binomial", 
                                Stat_Procedure = "Residual")

Stat_Comb_RF_C_IS=list(Fun_Stat = RF, model_type="classification", 
                       Stat_Procedure = "Direct")  

Stat_Comb_RF_C_Distill_RF=list(Fun_Stat = RF, model_type="classification", 
                               Fun_FitNull = RF_Null, null_model_type="classification", 
                               Stat_Procedure = "Distill") 

List_Stat_Comb_Logistic=c(
  'Stat_Comb_Logistic_Dev',    ## GLM-DEV     
  'Stat_Comb_RF_C_IS'          ## RF 
  
) 

List_Stat_Comb_Logistic_HD=c(
  "Stat_Comb_Logistic_Distill_Dev" ,   ## GLM-L1-D
  'Stat_Comb_Sel_Logistic_SST'         ## GLM-L1-R-SST
  , 'Stat_Comb_RF_C_IS'                ## RF
  
) 

# Non-GLM classification

List_Stat_Comb_Binary=c(
  'Stat_Comb_RF_C_IS',                ## RF
  'Stat_Comb_RF_C_Distill_RF',        ## RF-D
  'Stat_Comb_RF_R_Res_RF'             ## RF-RR
)

List_Stat_Comb_Binary_HD=c(  
  'Stat_Comb_RF_C_IS',                ## RF
  'Stat_Comb_RF_C_Distill_RF',        ## RF-D
  'Stat_Comb_RF_R_Res_RF'            ## RF-RR
) 

GetStatComb=function(setting){
  return(switch(setting, 
                'l_l_w'=List_Stat_Comb_Linear, 
                'l_l_m'=List_Stat_Comb_Nonlinear, 
                'l_gl_w'=List_Stat_Comb_Logistic, 
                'l_gl_m'=List_Stat_Comb_Binary, 
                'h_l_w'=List_Stat_Comb_Linear_HD, 
                'h_l_m'=List_Stat_Comb_Nonlinear_HD, 
                'h_gl_w'=List_Stat_Comb_Logistic_HD, 
                'h_gl_m'=List_Stat_Comb_Binary_HD
  ))
}

Run_all_test=function(x,y,setT,CSS_sample,setting, List_Stat_Comb=NULL, IncludeBaseline=T,...){
  
  cat('CSS testing begin;')
  CSS_test_result=list()
  
  
  
  pvalues=c()
  use.time=c()
  for(name_stat_comb in List_Stat_Comb){
    
    stat_comb=get(name_stat_comb)
    cat(name_stat_comb)
    start.time <- Sys.time()
    new.pvalue=CSSCRT(y, x, setT, CSS_samples=CSS_sample, 
                      Fun_Stat = stat_comb$Fun_Stat, model_type=stat_comb$model_type, 
                      Fun_FitNull=stat_comb$Fun_FitNull, null_model_type=stat_comb$null_model_type,
                      Stat_Procedure = stat_comb$Stat_Procedure,...)
    
    pvalues=c(pvalues,new.pvalue)
    
    print(new.pvalue)
    
    time.taken <- Sys.time() - start.time; cat('-Time used: ', time.taken,'s;')
    use.time=c(use.time, time.taken)
    cat('>')
  }
  names(pvalues)=names(use.time)=List_Stat_Comb
  
  
  cat('CSS testing end;\n')
  
  ## Setting the baseline methods for each model
  
  # Baseline methods
  
  if(IncludeBaseline==T){  # avoid running at this point
    baseline.result=Run_baseline_test(x,y,setT,setting)
    
    pvalues=c(pvalues,   baseline.result$pvalues)
    use.time=c(use.time, baseline.result$timings)
  }
  
  
  return(list(pvalues=pvalues, 
              use.time=use.time))
  
}

run_and_time <- function(expr, x, y, setT) {
  start_time <- Sys.time()
  result <- eval(expr)
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  return(list(result = result, time = elapsed_time))
}



Run_baseline_test=function(x,y,setT ,setting){
  
  
  # Store the function calls and their associated model types
  method_calls <- list(
    'l_l_' = list(
      ANOVA_GLM = quote(ANOVA_GLM(x, y, setT, model_type = 'gaussian')),
      dCRT_Lasso = quote(dCRT(x, y, setT, model_type = "Gaussian_lasso")),
      dCRT_RF = quote(dCRT(x, y, setT, model_type = "RF"))
    ),
    'l_gl' = list(
      ANOVA_GLM = quote(ANOVA_GLM(x, y, setT, model_type = 'binomial')),
      dCRT_Lasso = quote(dCRT(x, y, setT, model_type = "Binomial_lasso")),
      dCRT_RF = quote(dCRT(x, y, setT, model_type = "RF"))
    ),
    'h_l_' = list(
      DebiasLasso = quote(DL_hdi(x, y, setT, model_type = "gaussian")),
      CGM = quote(CGM(x, y, setT, model_type = "linear")),
      dCRT_Lasso = quote(dCRT(x, y, setT, model_type = "Gaussian_lasso")),
      dCRT_RF = quote(dCRT(x, y, setT, model_type = "RF"))
    ),
    'h_gl' = list(
      CGM = quote(CGM(x, y, setT, model_type = "logistic")),
      dCRT_Lasso = quote(dCRT(x, y, setT, model_type = "Binomial_lasso")),
      dCRT_RF = quote(dCRT(x, y, setT, model_type = "RF"))
    ), 
    'a_st' = list( #application, stock
      DebiasLasso = quote(DL_hdi(x, y, setT, model_type = "gaussian")),
      dCRT_Lasso = quote(dCRT(x, y, setT, model_type = "Gaussian_lasso")),
      dCRT_RF = quote(dCRT(x, y, setT, model_type = "RF"))
    ),
    'a_bc'= list(
      ChiSq = quote(Chisq_Multi(x, y, setT)),
      dCRT_Lasso = quote(dCRT(x, as.numeric(y), setT, model_type = "Gaussian_lasso")),
      dCRT_RF = quote(dCRT(x, as.numeric(y), setT, model_type = "RF"))
    )
    
  )
  
  # Identify the correct method calls based on the 'setting' variable
  selected_methods <- method_calls[[substr(setting, 1, 4)]]
  
  # Apply the function calls and record the results and timings
  results <- lapply(selected_methods, run_and_time, x, y, setT)
  
  # Extract pvalues and timings
  pvalues <- sapply(results, `[[`, "result")
  timings <- sapply(results, `[[`, "time")
  
  # Store results and timings in a named list
  baseline_results <- list(pvalues = pvalues, timings = timings)
}

### Conducts a conditional randomization test based on exchangeable sampling

CSSCRT <- function(Y, X, setT, CSS_samples = NULL, graph=NULL,
                   Fun_Stat, 
                   Stat_Procedure=c("Direct","Distill", "SelDist",
                                    "Residual","SelRes"),
                   Fun_FitNull=NULL, k_sel=NULL,
                   type='One-sided',
                   randomize=T,M=100,L=1, ...) {
  
  # Generate exchangeable samples if not provided
  if (is.null(CSS_samples)) {
    CSS_samples <- exchangeable_sampling(X, graph,I = setT, M=M, L=L,bReduce=T)
  }else{
    M=length(CSS_samples$X_copies)
  }
  
  X_copies <- CSS_samples$X_copies
  X_hub <- CSS_samples$X_hub
  
  nT = length(setT)
  new_setT=1:nT
  
  if(Stat_Procedure=="Direct"){
    T0 <- Fun_Stat(X,Y,setT,...)
    X_tilde=X
    Ttilde <- lapply(1:M, function(j){
      X_tilde[,setT]=X_copies[[j]];
      return(Fun_Stat(X_tilde,Y,setT,...))
    }
    )
  }else if(Stat_Procedure=="Residual"){
    Model_FitNull=Fun_FitNull(X[,-setT,drop=F],Y, k=-1,...)
    Residual_Null= as.numeric(Y) - as.numeric(Model_FitNull$fitted.value)
    
    T0 <- Fun_Stat(X[,setT,drop=F], Residual_Null,...)
    
    Ttilde <- lapply(1:M, function(j){
      Fun_Stat(X_copies[[j]],Residual_Null,...)
    })
    
  }else if(Stat_Procedure=="SelRes"){
    ## Specify the number of importance
    if (is.null(k_sel)){
      p=ncol(X)
      k_sel = min(as.integer(2 * log(p)), p - nT ) 
    }else{
      k_sel = 0  # Let the function itself choose
    }
    
    Model_FitNull=Fun_FitNull(X[,-setT,drop=F],Y, k_sel,...)
    Residual_Null= as.numeric(Y) - as.numeric(Model_FitNull$fitted.value)
    setS=Model_FitNull$imp.var
    k_sel=length(setS)
    
    X_TS=X[,c(setT,setS),drop=F]
    T0 <- Fun_Stat(X_ST,Residual_Null, setT=new_setT,...)
    Ttilde <- lapply(1:M, function(j){
      X_ST[,new_setT]=X_copies[[j]]
      Fun_Stat(X_ST,Residual_Null, setT=new_setT,...)
    })
    
  }else if(Stat_Procedure=="Distill"){
    
    Model_FitNull=Fun_FitNull(X[,-setT,drop=F],Y, k=-1,...)
    
    
    T0 <- Fun_Stat( cbind( X[,setT,drop=F], Model_FitNull$fitted.value),   
                    setT=new_setT, Y,...)
    
    
    Ttilde <- lapply(1:M, function(j){
      Fun_Stat(cbind(X_copies[[j]], Model_FitNull$fitted.value),   
               setT=new_setT, Y,...)
      
    })
    
  }else if(Stat_Procedure=="SelDist"){
    ## Specify the number of importance
    if (is.null(k_sel)){
      p=ncol(X)
      k_sel = min(as.integer(2 * log(p)), p - nT ) 
    }else{
      k_sel = 0  # Let the function itself choose
    }
    
    Model_FitNull=Fun_FitNull(X[,-setT,drop=F],Y, k_sel,...)
    
    setS=Model_FitNull$imp.var
    k_sel=length(setS)
    
    
    X_distill = cbind(X[,c(setT,setS),drop=F], Model_FitNull$fitted.value)
    
    
    T0 <- Fun_Stat(X_distill, Y, setT=new_setT,...)
    Ttilde <- lapply(1:M, function(j){
      X_distill[, new_setT]=X_copies[[j]]
      Fun_Stat(X_distill,Y, setT=new_setT,...)
    })
    
  }
  
  return(ComputePValue(T0, Ttilde, type=type, randomize = randomize))
}

### Conditional Randomization Test

# Function to simulate x
simulate_data <- function(Sigma, n) {
  X <- mvrnorm(n, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
  return(X)
}

# Function to simulate y
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

# Algorithm 4

CRT =function(popul_param,
              CSS_param=list(M=100,L=3),
              epo, 
              CSS_sample=NULL,  dir_prefix=''
){
  
  p <- popul_param$P
  n <- popul_param$N
  setT <- popul_param$setT
  
  setting <- popul_param$setting
  modelparam <- popul_param$modelParam
  
  #print(paste(n,p,setting,sep=","))
  path_prefix=paste0(dir_prefix,'DataStorage-',format(Sys.time(), "%b-%e-%Y"),'/')
  #     # Create the folder if the folder not exist
  if (!file.exists(path_prefix)) {
    dir.create(path_prefix)
  }
  
  file_prefix=paste0(c(rbind(names(popul_param)[1:3],unlist(popul_param[1:3]))),collapse ="_")
  
  file_prefix=paste0(path_prefix,file_prefix,collapse ="_")
  
  true_graph=graph_band(p,modelparam$K,modelparam$s)

  set.seed(epo) 
  perm.var=sample(p,p)
  sigma=true_graph$sigma[perm.var,perm.var]
  inv.sigma=true_graph$inv.sigma[perm.var,perm.var]
  
  diag_s=sqrt(diag(sigma))
  sigma=sweep(sweep(sigma,1,1/diag_s,"*"),2,1/diag_s, "*")
  inv.sigma=sweep(sweep(inv.sigma,1,diag_s,"*"),2,diag_s, "*")
  
  G = prec_to_adj(inv.sigma)
  
  if(min(eigen(inv.sigma)$values)<0){
    print('Model Error!')
    return(NULL)
  }
  if( 2*max(colSums(G)) + 3 > n){
    print('True graph too large!')
    return(NULL)
  }
  # Simulate the data
  x <- simulate_data(sigma, n)
  # saveRDS(x, file = paste0(file_prefix,"_n=",n,"_p=",p,"_x_epo=",epo,".rds") )
  beta=random_beta(p)
  noise_unif=runif(n)

  ### CSS sampling
  if(is.null(CSS_sample)){
    
    start.time <- Sys.time()
    #cat('CSS begin;')
    CSS_sample=exchangeable_sampling(x,G, I = setT, bReduced = T,
                                     M=CSS_param$M, L=CSS_param$L)
    time.taken <- Sys.time() - start.time
    #cat('Time used: ', time.taken,'s;')
    #cat('CSS end;\n')
    
    saveRDS(CSS_sample, file = paste0(file_prefix,"_n=",n,"_p=",p,"_CSS_epo=",epo,".rds") )
  }
  
  List_Stat_Comb=GetStatComb(setting)
  
  test_results=list()
  
  for(ind_theta in 1:length(modelparam$grid_signal)){
    
    theta=modelparam$grid_signal[ind_theta]
    #print(theta)
    
    y=simulate_y(x,setting,setT,theta,beta,
                 noise_unif=noise_unif)
  
    ### Inference begins
    #calculate statistics
    new_test_result=Run_all_test(x,y,setT,CSS_sample,setting, List_Stat_Comb)
    
    test_results[[ind_theta]]=list(
      theta=theta,
      pvalues=new_test_result$pvalues,
      use.time=new_test_result$use.time
    )
    #print(unlist(test_results[[ind_theta]]))
  }
  
  return(test_results)
}     
