
wd_temp=getwd()
setwd('../../CRT')
source('../GGMTesting.R')
source('../CSS-GGM.R')
source('CRT-Stat.R')
source('CRT-Computing.R')
source('Run-Testing.R')
setwd(wd_temp)


para_run_one=function(X,Y,setT,G1){

  
  List_Stat_Comb = c(
    "Stat_Comb_Lasso_Res_Linear_Dev", 
    "Stat_Comb_RF_R_Res_RF"
  )
  
  setting='a_st'   # the setting is application: stock
  
CSS_param=list(M=400,L=3)
CSS_sample=exchangeable_sampling(X,G1, I = setT, bReduced = T,
                                 M=CSS_param$M, L=CSS_param$L)

CI_test_result=Run_all_test(X,Y,setT,CSS_sample,setting = setting, List_Stat_Comb = List_Stat_Comb,
                            IncludeBaseline=T )

print(CI_test_result)
}