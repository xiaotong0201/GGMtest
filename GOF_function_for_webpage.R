

# GGM_webpage -------------------------------------------------------------

GGM_webpage =function(popul_param,
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
  CSS_test_T0 = unlist(lapply(CSS_test_result, function(a)a$T0))
  CSS_test_Ttilde = map(CSS_test_result,.f = function(x){return(x$Ttilde)})
  names(CSS_test_pvalues)=List_CSS_stat
  names(CSS_test_T0)=List_CSS_stat
  names(CSS_test_Ttilde)=List_CSS_stat
  # Wrap results in a list
  final_result <- list(
    p_value = list(
      'VV_pvalue' = VV_pvalue,
      'DP_pvalue' = DP_pvalue,
      'CSS_test_pvalues' = CSS_test_pvalues),
    T0 = CSS_test_T0,
    Ttilde = CSS_test_Ttilde
  )
  return(final_result)
}     

# run GGM experiments for webpage
run_experiments_epo_webpage <- function(param_template, epo, M,L,List_CSS_stat=NULL,
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
        final_results[[lFR+1]] = GGM_webpage(param_part, 
                                             CSS_param=list(M=M,L=L),
                                             epo=epo, List_CSS_stat=List_CSS_stat,
                                             dir_prefix=dir_prefix)
        lFR=lFR+1
      }
    }
  }
  return(final_results)
}

GoF_test_webpage = function(X, X_copies, G0, List_CSS_stat = NULL){
  
  if(is.null(List_CSS_stat)==T){
    
    List_CSS_stat = c( "PRC_SS", "ERC_Z_SS", "F_sum", "glr_glasso")
    
  }
  
  CSS_test_result = list()
  i = 1
  while(i<= length(List_CSS_stat)){
    f=List_CSS_stat[i]
    
    if(substr(f,2,3)=="RC"){
      temp_fn=function(X,g)All_Stat_RC(X,g,residual_type=substr(f,1,2))
      bundle_results=CSSGoF(X, G0, X_copies, Fun_Stat = temp_fn)
      CSS_test_result[[i]] = list(
        pvalue=bundle_results$pvalue[f],
        T0=bundle_results$T0[f],
        Ttilde=bundle_results$Ttilde[f,]
      )
      
    }else if(startsWith(f,"F_")){
      bundle_results=CSSGoF(X,G0, X_copies, Fun_Stat = All_Stat_F)
      
      CSS_test_result[[i]] = list(
        pvalue=bundle_results$pvalue[f],
        T0=bundle_results$T0[f],
        Ttilde=bundle_results$Ttilde[f,]
      )
      
    }else{
      CSS_test_result[[i]] = CSSGoF(X,G0, X_copies, Fun_Stat = get(f))
    }
    
    i=i+1
  }
  
  CSS_test_pvalues = unlist(lapply(CSS_test_result, function(a)a$pvalue))
  CSS_test_T0 = unlist(lapply(CSS_test_result, function(a)a$T0))
  CSS_test_Ttilde = map(CSS_test_result,.f = function(x){return(x$Ttilde)})
  names(CSS_test_pvalues)=List_CSS_stat
  names(CSS_test_T0)=List_CSS_stat
  names(CSS_test_Ttilde)=List_CSS_stat
  
  GoF_test_result = list(
    p_values = CSS_test_pvalues,
    T0s = CSS_test_T0,
    Ttildes = CSS_test_Ttilde
  )
  
  return(GoF_test_result)
}