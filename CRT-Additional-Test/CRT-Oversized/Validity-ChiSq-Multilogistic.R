p <- 10
n <- 50  #500 or 50
setT <- 1:4

modelparam <- 'l-gl-w'
library(MASS)
wd=getwd()
setwd('../../CRT-Oversized')
source('CRT-Simulate.R')
setwd(wd)

true_graph=graph_band(p, 5,  0.2 )
pv=c()

library(VGAM)

Chisq_multi=function(X,Y,setT){
  
  
  reduced = vglm(Y ~ X[, -setT,drop=F],family = multinomial)
  full = vglm(Y ~ X, family = multinomial)
  
  
  lrt=lrtest(full, reduced)
  
  return(lrt@Body[2,5])
}

beta1=random_beta(p)
beta2=random_beta(p)
beta3=random_beta(p)

for(epo in 1:800){
  set.seed(epo) 
  perm.var=sample(p,p)
  sigma=true_graph$sigma[perm.var,perm.var]
  inv.sigma=true_graph$inv.sigma[perm.var,perm.var]
  
  diag_s=sqrt(diag(sigma))
  sigma=sweep(sweep(sigma,1,1/diag_s,"*"),2,1/diag_s, "*")
  inv.sigma=sweep(sweep(inv.sigma,1,diag_s,"*"),2,diag_s, "*")
  
  G = prec_to_adj(inv.sigma)
  
  # Simulate the data
  x <- simulate_data(sigma, n)
 
  noise_unif=runif(n)
  
  fit_value_null=exp(x[,-setT]%*%cbind(beta1[-setT],beta2[-setT],beta3[-setT] ))
  
  
  Y_sim=t(apply(fit_value_null, 1, function(v)rmultinom(1,1,v/sum(v))))
  Y_sim_factor=apply(Y_sim==1,1,which)
  
  
  pv[epo]=Chisq_multi(x, Y_sim, setT)
  
}

hist(pv)
