
# Run_CIT_testing_webpage -----------------------------------------------------

# run_all_test function for webpage

Run_all_test_webpage=function(x,y,setT,CSS_sample,setting, List_Stat_Comb=NULL, IncludeBaseline=T){
  
  #cat('CSS testing begin;')
  CSS_test_result=list()
  
  pvalues=c()
  use.time=c()
  T0s = list()
  Ttildes = list()
  for(name_stat_comb in List_Stat_Comb){
    
    stat_comb=get(name_stat_comb)
    #cat(name_stat_comb)
    start.time <- Sys.time()
    new.pvalue=CSSCRT_webpage(y, x, setT, CSS_samples=CSS_sample, 
                      Fun_Stat = stat_comb$Fun_Stat, model_type=stat_comb$model_type, 
                      Fun_FitNull=stat_comb$Fun_FitNull, null_model_type=stat_comb$null_model_type,
                      Stat_Procedure = stat_comb$Stat_Procedure)$pvalue
    
    pvalues=c(pvalues,new.pvalue)
    
    #print(new.pvalue)
    
    T0 = CSSCRT_webpage(y, x, setT, CSS_samples=CSS_sample, 
                Fun_Stat = stat_comb$Fun_Stat, model_type=stat_comb$model_type, 
                Fun_FitNull=stat_comb$Fun_FitNull, null_model_type=stat_comb$null_model_type,
                Stat_Procedure = stat_comb$Stat_Procedure)$T0
    
    Ttilde = unlist(CSSCRT_webpage(y, x, setT, CSS_samples=CSS_sample, 
                           Fun_Stat = stat_comb$Fun_Stat, model_type=stat_comb$model_type, 
                           Fun_FitNull=stat_comb$Fun_FitNull, null_model_type=stat_comb$null_model_type,
                           Stat_Procedure = stat_comb$Stat_Procedure)$Ttilde)
    name <- name_stat_comb
    T0s[[name]]=T0
    Ttildes[[name]]=Ttilde
    
    time.taken <- Sys.time() - start.time
    #cat('-Time used: ', time.taken,'s;')
    use.time=c(use.time, time.taken)
    #cat('>')
  }
  names(pvalues)=names(use.time)=List_Stat_Comb
  
  
  #cat('CSS testing end;\n')
  if(IncludeBaseline==T){  # avoid running at this point
    baseline.result=Run_baseline_test(x,y,setT,setting)
    
    pvalues=c(pvalues,   baseline.result$pvalues)
    use.time=c(use.time, baseline.result$timings)
  }
  
  
  return(list(pvalues=pvalues,
              T0s = T0s,
              Ttildes = Ttildes,
              use.time=use.time))
  
}


# CSSCRT_webpage ---------------------------------------------------------


CSSCRT_webpage <- function(Y, X, setT, CSS_samples = NULL, graph=NULL,
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
  
  
  return(list(
    pvalue = ComputePValue(T0, Ttilde, type=type, randomize = randomize),
    T0 = T0,
    Ttilde = Ttilde))
}


# CRT_webpage -------------------------------------------------------------

CRT_webpage =function(popul_param,
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
    new_test_result=Run_all_test_webpage(x,y,setT,CSS_sample,setting, List_Stat_Comb)
    
    
    test_results[[ind_theta]]=list(
      theta=theta,
      pvalues=new_test_result$pvalues,
      T0s = new_test_result$T0s,
      Ttildes = new_test_result$Ttildes,
      use.time=new_test_result$use.time
    )
    #print(unlist(test_results[[ind_theta]]))
  }
  
  return(test_results)
}     

GRT_test_webpage= function(X, CSS_sample, setting, setT, beta, niose_unif, grid_signal, List_Stat_Comb,IncludeBaseline = T){
  
  test_results=list()
  
  for(ind_theta in 1:length(grid_signal)){
    
    theta=grid_signal[ind_theta]
    
    y=simulate_y(X,setting,setT,theta,beta,
                 noise_unif=noise_unif)
    
    
    ### Inference begins
    #calculate statistics
    
    
    pvalues=c()
    use.time=c()
    T0s = list()
    Ttildes = list()
    for(name_stat_comb in List_Stat_Comb){
      
      stat_comb=get(name_stat_comb)
      
      new.pvalue=CSSCRT_webpage(y, X, setT, CSS_samples=CSS_sample, 
                                Fun_Stat = stat_comb$Fun_Stat, model_type=stat_comb$model_type, 
                                Fun_FitNull=stat_comb$Fun_FitNull, null_model_type=stat_comb$null_model_type,
                                Stat_Procedure = stat_comb$Stat_Procedure)$pvalue
      
      pvalues=c(pvalues,new.pvalue)
      
      T0 = CSSCRT_webpage(y, X, setT, CSS_samples=CSS_sample, 
                          Fun_Stat = stat_comb$Fun_Stat, model_type=stat_comb$model_type, 
                          Fun_FitNull=stat_comb$Fun_FitNull, null_model_type=stat_comb$null_model_type,
                          Stat_Procedure = stat_comb$Stat_Procedure)$T0
      
      Ttilde = unlist(CSSCRT_webpage(y, X, setT, CSS_samples=CSS_sample, 
                                     Fun_Stat = stat_comb$Fun_Stat, model_type=stat_comb$model_type, 
                                     Fun_FitNull=stat_comb$Fun_FitNull, null_model_type=stat_comb$null_model_type,
                                     Stat_Procedure = stat_comb$Stat_Procedure)$Ttilde)
      name <- name_stat_comb
      T0s[[name]]=T0
      Ttildes[[name]]=Ttilde
    }
    names(pvalues)=List_Stat_Comb
    if(IncludeBaseline==T){  
      baseline.result=Run_baseline_test(X,y,setT,setting)
      
      pvalues=c(pvalues,   baseline.result$pvalues)
      
    }
    new_test_result= list(pvalues=pvalues,
                          T0s = T0s,
                          Ttildes = Ttildes)
    
    name = paste("theta=",theta)
    test_results[[name]]=list(
      theta=theta,
      pvalues=new_test_result$pvalues,
      T0s = new_test_result$T0s,
      Ttildes = new_test_result$Ttildes
    )
  }
  return(test_results)
  
}
