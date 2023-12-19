
source('CRT-Stat.R')

## Setting the test statistics for each model

# Linear regression----

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


# Non-Linear regression-------


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


# Logistic regression  ---------

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

# Non-GLM classification  ----

List_Stat_Comb_Binary=c(
  'Stat_Comb_RF_C_IS',                ## RF
  'Stat_Comb_RF_C_Distill_RF',        ## RF-D
  'Stat_Comb_RF_R_Res_RF',            ## RF-RR
  )

List_Stat_Comb_Binary_HD=c(  
  'Stat_Comb_RF_C_IS',                ## RF
  'Stat_Comb_RF_C_Distill_RF',        ## RF-D
  'Stat_Comb_RF_R_Res_RF',            ## RF-RR
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

# Baseline methods ----

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


