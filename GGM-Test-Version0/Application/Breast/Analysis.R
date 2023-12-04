# Y part 
# Load the CSV data
table_path = 'GSE2034-PI.csv'
PI = read.csv(table_path)

# Display the first few rows of the loaded data
head(PI)

# subclass
ER=factor(PI$ER.Status,levels=c("ER+","ER-"))

# Survival time
ST=PI$time.to.relapse.or.last.follow.up..months.
NotCensor=PI$relapse..1.True.

# Create a binary variable Y.brain based on conditions
Y.brain = rep(0, nrow(PI))
Y.brain[which(PI[, 4] <= 48)] = 1
Y.brain[which(PI[, 4] <= 48 & PI[, 7])] = 2
Y.brain = as.factor(Y.brain)

# Display a summary of Y.brain
summary(Y.brain)
Y=Y.brain

# X part 
# Load gene data
Data = read.csv('GSE2034_series_matrix.txt', skip = 26, sep = '\t')

GeneData = t(Data[-(1:28),])
GeneName = GeneData[1,]
GeneData = apply(GeneData[-1,],MARGIN = c(1, 2), FUN = as.numeric)

ImpVarName = read.csv('ImpVar.csv',header = F)$V1
ImpVarInd =  which(GeneName %in% ImpVarName)  #match(ImpVarName,GeneName) #
X_gene = GeneData[
  match(PI$GEO.asscession.number,Data[1,-1])
  ,ImpVarInd]


# Display the dimensions of X
dim(X_gene)

# Create a histogram for the first column of X
hist(X_gene[,2])


library(bestNormalize)
X_norm=apply(X_gene, 2, function(v)orderNorm(v)$x.t)
order_shapiro_p_values <-
  apply(X_norm, 2, function(column)
    shapiro.test(column)$p.value)
sort(order_shapiro_p_values, index=T)

#############################
### check if the mariginal is significant

X=X_norm

for(j in 1:ncol(X)){
  linearmodel=glm((Y!=0)~X[,j],family=binomial)
  print(summary(linearmodel)$coef[2,4])
}

for(j in 1:ncol(X)){
  cox_model=coxph(Surv(time = ST,event = NotCensor) ~ X[,j])
  print(summary(cox_model)$coef[1,5])
}

##############################

library(glmnet)
# Goal: Define T  
all.rel.set=1:76

# May use ER
sum(ER=='ER+')

sub_X=X[ER=='ER+', 1:60 ]
sub_Y=Y[ER=='ER+']


# Method 1, use all stocks to run a Lasso selection
set.seed(2023)
lasso.cv=cv.glmnet(sub_X,sub_Y,family="multinomial",nfolds = 5)
lasso=glmnet(sub_X ,sub_Y,family="multinomial")
betahat_list=coef(lasso,s=lasso.cv$lambda.min)
betahat_mat=do.call(cbind,betahat_list)
betahat_nonzero=which((rowSums(betahat_mat^2)[-1])!=0)
strong.set=intersect(all.rel.set, betahat_nonzero)
weak.set=setdiff(all.rel.set,betahat_nonzero)


## Meaning: whether should use the best CV lambda or should 
##  include more predictors
setT_multin_weak=weak.set


plot(lasso)
betahat_ose_list=coef(lasso,s=lasso.cv$lambda.1se)
betahat_ose_mat=do.call(cbind,betahat_ose_list)
betahat_ose_nonzero=which((rowSums(betahat_ose_mat^2)[-1])!=0)

## Meaning: whether should use the best CV lambda or the 1se lambda
setT_multin_1se=setdiff(strong.set,betahat_ose_nonzero) 


GeneName[ImpVarInd[betahat_ose_nonzero]]
GeneName[ImpVarInd[setT_multin_1se]]

######### CI testing ################


wd_temp = getwd()
setwd('../../CRT')
source('../CSS-GGM.R')
source('../GGMTesting.R')
source('CRT-Stat.R')
source('CRT-Computing.R')
source('Run-Testing.R')
setwd(wd_temp)



############
mat_X=sub_X[, strong.set]
p=ncol(mat_X)
setT= match(setT_multin_1se , strong.set)


Chisq_Multi(mat_X, sub_Y, setT)


G1=matrix(1,p,p)
diag(G1)=0




exp_seed=2023
set.seed(exp_seed)

  CSS_param = list(M = 1000, L = 3)
  
  CSS_sample = exchangeable_sampling(
    mat_X,
    G1,
    I = setT,
    bReduced = T,
    M = CSS_param$M,
    L = CSS_param$L
  )
  
  ## random forest ok
  RF_result = Run_all_test(
    mat_X,
    sub_Y, 
    setT,
    CSS_sample,
    setting = 'a_bc',
    List_Stat_Comb = "Stat_Comb_RF_C_IS" ,
    IncludeBaseline = T, 
  )
  print(RF_result)
  
  save(mat_X, sub_Y, setT, G1,CSS_sample, RF_result, file='ERplus-CI-testing-seed2023.RData' )
  
