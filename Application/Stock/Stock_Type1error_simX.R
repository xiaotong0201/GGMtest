library(parallel)
### type i error check
source('~/Library/CloudStorage/GoogleDrive-statdongminghuang@gmail.com/My Drive/Codes in R/Template_foreach_evendistributed.R')


load("DataForExp.RData")

prefunc = function() {
  setwd('../../CRT')
  source('../GGMTesting.R')
  source('../CSS-GGM.R')
  source('CRT-Stat.R')
  source('CRT-Computing.R')
  source('Run-Testing.R')
}


############# baseline

baseline_method_calls = list(
  DebiasLasso = quote(DL_hdi(x, y, setT, model_type = "gaussian")),
  dCRT_Lasso = quote(dCRT(x, y, setT, model_type = "Gaussian_lasso")),
  dCRT_RF = quote(dCRT(x, y, setT, model_type = "RF"))
)



null_result = list()
null_func = function(i) {
  set.seed(i)
  
  # Generate simulated data
  sim_sample = exchangeable_sampling(
    mat_X,
    G1, 
    bReduced = F,
    M = 1,
    L = 20
  )
  sim_X=sim_sample$X_copies[[1]] 
  Y_sim =coef_null[1]+sim_X[,-setT]%*% c(coef_null[-1]) + err_sd_null * rnorm(nrow(mat_X))
  
  
  CSS_param = list(M = 400, L = 3)
  CSS_sample = exchangeable_sampling(
    sim_X,
    G1,
    I = setT,
    bReduced = T,
    M = CSS_param$M,
    L = CSS_param$L
  )
  
  CI_test_result = Run_all_test(
    sim_X,
    Y_sim,
    setT,
    CSS_sample,
    setting = 'a_st',
    List_Stat_Comb = List_Stat_Comb,
    IncludeBaseline = T
  )
  return(CI_test_result)
}

List_Stat_Comb = c(
  "Stat_Comb_Lasso_Res_Linear_Dev",
  "Stat_Comb_RF_R_Res_RF"
)



setT=setT2

setT=setT1


setV=setdiff(1:ncol(mat_X), setT)
cv_lasso_null=cv.glmnet(mat_X[, setV], Y)
fit_value_null=predict(cv_lasso_null,newx = mat_X[, setV], type='response', s=cv_lasso_null$lambda.min)
err_sd_null=sd(Y-fit_value_null)
coef_null = predict(cv_lasso_null,type='coef', s=cv_lasso_null$lambda.min)


null_result = run_parallel(
  replication_start = 1,
  replication_end = 400,
  num_cores = 8,
  list_export = c(
    'err_sd_null',
    'coef_null',
    'mat_X',
    'G1',
    'run_and_time',
    'baseline_method_calls',
    'setT', 'List_Stat_Comb'
  ),
  func_name = 'null_func',
  prefunc_name = 'prefunc'
)

if(sum(setT!=setT2)==0){
  save(null_result, file='Type-I_Strong4_DIA_2023-11-12.RData')
}else if(sum(setT!=setT1)==0){
  save(null_result, file='Type-I_Weak5_DIA_2023-11-12.RData')
}


