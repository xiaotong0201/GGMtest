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



param_VV_template <- list(model = "Verzelen",
                          N = c(10,15,30),
                          P = c(15),
                          modelParam = list(
                            Q = c(0.1,0.4,1),
                            E = c(0.1,0.15)
                          )
)


run_experiments_epo <- function(param_template, epo, M,L,List_CSS_stat=NULL,
                                dir_prefix='',file) {
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
        
        save(final_results,param_template,
             file=file)
        
      }
    }
  }
  
  return(final_results)
}



args=commandArgs(trailingOnly = TRUE)
ind.experiment=as.numeric(args[1])
seq.id = as.numeric(args[2])



FlagLocalTest=F

# small test
if(FlagLocalTest){
  
  L <- 2
  M <- 20 
  
  
  param_test_template <- list(model = "Verzelen",
                            N = c(10),
                            P = c(15),
                            modelParam = list(
                              Q = c(0.1),
                              E = c(0.1)
                            )
  )
  
  param_template=param_test_template
  param_template=param_VV_template
  ind.experiment=4
  seq.id = 1
}else{
  
  param_template=param_VV_template
  
  
  L=3
  M=100
}

PlatForm=(Sys.info())["sysname"]
if(FlagLocalTest&(PlatForm!="Windows")){Rprof("my_profile.out")}


path_prefix=paste0('../../../Programming/LocalExperiments/GGM/GoF/Experiments-',ind.experiment,'/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

for( epo in (seq.id-1)*16+1:16){

save(PlatForm,file=
       paste0(path_prefix,epo,'.RData')
     )

###############################################################################


output = run_experiments_epo(param_template,epo, M,L,dir_prefix=path_prefix,
                             List_CSS_stat=c(
                               "PRC_SS", "PRC_SA", "PRC_Z_SS", "PRC_Z_SA",
                               "ERC_SS", "ERC_SA", "ERC_Z_SS", "ERC_Z_SA",
                               "F_max", "F_sum","F_logsum", "F_Z_SA","F_Z_SS"
                             ),paste0(path_prefix,epo,'.RData'))
print(output)

save(output,param_template,
     file=paste0(path_prefix,epo,'.RData'))

###############################################################################


if(FlagLocalTest&(PlatForm!="Windows")){
  Rprof(NULL)
  summary_data <- summaryRprof("my_profile.out")
  print(summary_data)
}


}