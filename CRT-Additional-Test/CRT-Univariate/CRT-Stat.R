library(glmnet)
library(rags2ridges)
library(glasso)



## For direct computations ####

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





## For two-step computation ####

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


###################### Baseline Functions #################

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
  dCRT_res = d_CRT(X,Y,FDR = 1, model = model_type,CDR='No',candidate_set=setT)
  return(bonfTest(dCRT_res$p_values))   # a bug fixed on 2023-11-07 
}



