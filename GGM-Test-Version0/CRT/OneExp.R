## Graphical conditional randomization test (G-CRT)
## Output the p-values for each statistic at different values of theta.

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
  
  
  
  print(paste(n,p,setting,sep=","))
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
    cat('CSS begin;')
    CSS_sample=exchangeable_sampling(x,G, I = setT, bReduced = T,
                                     M=CSS_param$M, L=CSS_param$L)
    time.taken <- Sys.time() - start.time; cat('Time used: ', time.taken,'s;')
    cat('CSS end;\n')
    
    
    saveRDS(CSS_sample, file = paste0(file_prefix,"_n=",n,"_p=",p,"_CSS_epo=",epo,".rds") )
    }
  
  List_Stat_Comb=GetStatComb(setting)
  
  
  test_results=list()
  
  for(ind_theta in 1:length(modelparam$grid_signal)){
    
    theta=modelparam$grid_signal[ind_theta]
    print(theta)
    
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
    print(unlist(test_results[[ind_theta]]))
  }
  
  return(test_results)
}     




