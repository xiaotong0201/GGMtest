rm(list=ls())

FlagLocalTest=F

library(Matrix)
library(MASS)
library(mvtnorm)
source('../CSS-GGM.R')
source('../GGMTesting.R')
source('CRT-Stat.R')
source('CRT-Simulate.R')
source('CRT-Computing.R')
source('Run-Testing.R')



source('OneExp.R')

## Low dimension----
popul_param_llw <- list(
  setting = "l_l_w",
  N = c(50),
  P = c(20),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= seq(0,2.5,by=0.25)
  )
)

popul_param_hlw <- list(
  setting = "h_l_w",
  N = c(80),
  P = c(120),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= seq(0, 1.5, by=0.25)
  )
)

popul_param_lgw <- list(
  setting = "l_gl_w",
  N = c(50),
  P = c(20),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal = c(0,1,2,3, 4,6, 8) # c(0,1,2,3,5,7.5,10) #seq(0,3,by=0.5)
  )
)


popul_param_hgw <- list(
  setting = "h_gl_w",
  N = c(80),
  P = c(120),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= seq(0, 3.5, by=0.5) ### c(0,1,2,3,5,7.5,10)
  )
)



popul_param_llm <- list(
  setting = "l_l_m",
  N = c(50),
  P = c(20),
  setT = 1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= c(0,1,2,3,5,7) # c(0,1,2,3,5,7.5,10,12.5,15) #seq(0,10,by=1)
  )
)



popul_param_hlm <- list(
  setting = "h_l_m",
  N = c(80),
  P = c(120),
  setT = 1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= c(0,1,2,3,5,7) # c(0,1,2,3,5,7.5,10,12.5,15)
  )
)



popul_param_lgm <- list(
  setting = "l_gl_m",
  N = c(50),
  P = c(20),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= c(0,1,2,3,5,7)
  )
)

popul_param_hgm <- list(
  setting = "h_gl_m",
  N = c(80),
  P = c(120),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= c(0,1,2,3,5,7)
  )
)




args=commandArgs(trailingOnly = TRUE)
ind.experiment=as.numeric(args[1])
epo = as.numeric(args[2])




# small test
if(FlagLocalTest){
  
  popul_param=list(
    setting = "h_l_w",
    N = c(80),
    P = c(120),
    setT=1, 
    modelParam = list(
      K = 6,
      s = 0.2,
      grid_signal= seq(0, 1.5, by=0.25)
    )
  )
  epo=24
  ind.experiment=5
  
  
}else{
  List_Settings=c('llw','lgw','llm','lgm',
                  'hlw','hgw','hlm','hgm')
  
  popul_param=get(paste0('popul_param_',List_Settings[ind.experiment]))
  
}

PlatForm=(Sys.info())["sysname"]
if(FlagLocalTest&(PlatForm!="Windows")){Rprof("my_profile.out")}


path_prefix=paste0('/Users/stahd/GitHub/Programming/LocalExperiments/GGM/CRT-TDistribution/Experiments-',ind.experiment,'/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

save(PlatForm,file=paste0(path_prefix,epo,'.RData'))

###############################################################################


output = CRT(popul_param,epo = epo)

save(output,popul_param,
     file=paste0(path_prefix,epo,'.RData'))

###############################################################################


if(FlagLocalTest&(PlatForm!="Windows")){
  Rprof(NULL)
  summary_data <- summaryRprof("my_profile.out")
  print(summary_data)
}


