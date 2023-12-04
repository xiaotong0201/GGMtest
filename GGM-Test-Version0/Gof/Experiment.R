rm(list=ls())

library(Matrix)
library(MASS)
library(mvtnorm)
library(glasso)
library(Rfast)
library(glmnet)
library(rags2ridges)
source('../CSS-GGM.R')
source('../GGMTesting.R')
source('GoF-Statistics.R')
source('GoF-Settings.R')
source('GoF-OneExp.R')



FlagLocalTest=F



param_Band_Null_template <- list(
  model = "Band_Band",
  N = c(20, 40, 80),
  P = c(20, 120),
  modelParam = list(
    K0 = 6,
    K = 6,
    s = 0.2 
  )
)  ### 50 minutes for one



param_Band_template <- list(
  model = "Band_Band",
  N = c(20, 40, 80),
  P = c(20, 120),
  modelParam = list(
    K0 = 1,
    K = c(2,  6),
    s = c(0.1, 0.15,0.2)
  )
)




param_Hub_template <- list(
  model = "Hub",
  N = c(20, 40, 80),
  P = c(20, 120),
  modelParam = list(hubsize=10,
                    Q = 0.7,
                    noise = c(0.5, 0.9, 1.4, 2))
)


param_ER_template <- list(model = "ER",
                          N = c(50,100),
                          P = c(40,120),
                          modelParam = list(
                            Q = c(0.2,0.4),
                            Q0 = 0.08,
                            s = c(0.01,0.02)
                          )
)


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



args=commandArgs(trailingOnly = TRUE)
ind.experiment=as.numeric(args[1])
epo = as.numeric(args[2])



# small test
if(FlagLocalTest){
  
  L <- 2
  M <- 20 
  (1-pbinom( (50-3)/2,120,0.08,lower.tail=F))^200
  pbinom( (40-3)/2,20,0.12,lower.tail=F)
  
  param_test_template <-  list(model = "ER",
                               N = c(50,100),
                               P = c(40,120),
                               modelParam = list(
                                 Q = c(0.2,0.4),
                                 Q0 = 0.08,
                                 s = c(0.01,0.02)
                               )
  )
  
  param_template=param_test_template
  ind.experiment=0
}else{
  
if(ind.experiment==0){
  param_template=param_Band_Null_template
}else if(ind.experiment==1){
    param_template=param_Band_template
  }else if(ind.experiment==2){
    param_template=param_Hub_template
  }else if(ind.experiment==3){
    param_template=param_ER_template
  }
  
  L=3
  M=100
}

PlatForm=(Sys.info())["sysname"]
if(FlagLocalTest&(PlatForm!="Windows")){Rprof("my_profile.out")}


path_prefix=paste0('../../../Programming/LocalExperiments/GGM/GoF/Experiments-',ind.experiment,'/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

save(PlatForm,file=paste0(path_prefix,epo,'.RData'))

###############################################################################


output = run_experiments_epo(param_template,epo, M,L,dir_prefix=path_prefix)
print(output)

save(output,param_template,
     file=paste0(path_prefix,epo,'.RData'))

###############################################################################


if(FlagLocalTest&(PlatForm!="Windows")){
  Rprof(NULL)
  summary_data <- summaryRprof("my_profile.out")
  print(summary_data)
}


