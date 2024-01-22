
### We summariz the setup for generating data and running GoF test of the simulation studies

source('Function.R')
library(Matrix)
library(MASS)
library(mvtnorm)
library(glasso)
library(Rfast)
library(glmnet)
library(rags2ridges)
library(dplyr)
library(reshape)
library(ggplot2)

##########################################################################
               ###### Type I error control ######
##########################################################################

# band graph (K,K_0,p,s,n)=(6,6,20, 0.2, 80) 

# Parameters
param_Band_template <- list(
  model = "Band_Band",
  N = 80,
  P = 20,
  modelParam = list(
    K0 = 6,
    K = 6,
    s = 0.2
  )
)
M = 100
L = 3

# Set the storage path
path_prefix=paste0('Experiments-type_I_error','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# Repeat experiments 1200 times and save results
for (epo in 1:1200){
  output = run_experiments_epo(param_Band_template,epo, M,L,List_CSS_stat =NULL ,dir_prefix=path_prefix)
  save(output,param_Band_template,file=paste0(path_prefix,epo,'.RData'))
}

# Load results
all_output=list()
for(epo in 1:1200){
  load(file=paste0('Experiments-type_I_error','/',epo,'.RData'))
  all_output[[epo]]=output
}

# calculate the rejection proportions and conduct the one-sample proportion test for the hypothesis that the rejection probability is no greater than alpha = 0.05.
# utilize four test statistic functions (PRC, ERC, F_Sigma, GLR-l1) and take M1P1 and Bonf as baseline methods.

summ=read_experiments_epo(param_Band_template,all_output)
#Rejection proportion
rej.table=extract_metric(metric_name='power')[c(1:8,14,17,22)]
colnames(rej.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_Sigma","GLR_l1")
#One-sided p-value of the one-sample proportion test
ospt.table=extract_metric(metric_name='ospt')[c(1:8,14,17,22)]
colnames(ospt.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_Sigma","GLR_l1")


# QQ plot comparing the uniform distribution and the p-values of CSS-GLR with GLR-l1.

# Create the QQ plot using ggplot2's stat_qq function
pvs=summ[[1]]$pvalue.table[,'glr_glasso']
qq_plot <- ggplot(data.frame(pvs), aes(sample = pvs)) +
  stat_qq(distribution = qunif, color='black', shape=1, size=0.1) + # Uniform distribution quantiles
  stat_qq_line(distribution = qunif, linetype = "dashed",color='red',linewidth=1) +
  theme_bw() +
  labs(
    x = "Theoretical Quantiles",
    y = "Sample P-Values"
    # ,
    # title='P-values of GLR-L1'
  ) +
  theme(
    text = element_text(size = 16), # Adjusts overall text size; you can change the value as needed
    axis.title.x = element_text(size = 18), # Adjusts x-axis label size
    axis.title.y = element_text(size = 18), # Adjusts y-axis label size
    plot.title = element_text(size = 20) # Adjusts plot title size
  )

# Print the plot
print(qq_plot)


##########################################################################
                  ###### Power comparison ######
##########################################################################


############### band graph (K=6,K_0=1,p=c(20,120),s=c(0.1,0.15,0.2),n=c(20,40,80)) ##########

param_Band_template <- list(
  model = "Band_Band",
  N = c(20,40,80),
  P = c(20,120),
  modelParam = list(
    K0 = 1,
    K = 6,
    s = c(0.1,0.15,0.2)
  )
)
M = 100
L = 3

path_prefix=paste0('Experiments-power-band','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# repeat experiment 400 times and save the result
for (epo in 1:400){
  output = run_experiments_epo(param_Band_template,epo, M,L,List_CSS_stat =NULL ,dir_prefix=path_prefix)
  save(output,param_Band_template,file=paste0(path_prefix,epo,'.RData'))
}

#load the result
all_output=list()
for(epo in 1:400){
  load(file=paste0('Experiments-power-band','/',epo,'.RData'))
  all_output[[epo]]=output
}

# Calculate the power and standard error
## utilize four test statistic functions (PRC, ERC, F_Sigma, GLR-l1) and take M1P1 and Bonf as baseline methods.

summ=read_experiments_epo(param_Band_template,all_output)
#Power
power.table=extract_metric(metric_name='power')[c(1:8,14,17,22)]
colnames(power.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_sum","GLR_lasso")
#Standard error
SE.table=extract_metric(metric_name='SE')[c(1:8,14,17,22)]
colnames(SE.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_sum","GLR_lasso")



############### hub graph (hubsize=10,Q=0.3,p=c(20,120),\xi=c(0.9,1.4,2),n=c(20,40,80)) ##########

param_Hub_template <- list(
  model = "Hub",
  N = c(20,40,80),
  P = c(20,120),
  modelParam = list(
    hubsize = 10,
    Q = 0.3,
    noise = c(0.9,1.4,2)
  )
)
M = 100
L = 3

path_prefix=paste0('Experiments-power-hub','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# repeat experiment 400 times and save the result
for (epo in 1:400){
  output = run_experiments_epo(param_Hub_template,epo, M,L,List_CSS_stat =NULL ,dir_prefix=path_prefix)
  save(output,param_Hub_template,file=paste0(path_prefix,epo,'.RData'))
}

#load the result
all_output=list()
for(epo in 1:400){
  load(file=paste0('Experiments-power-hub','/',epo,'.RData'))
  all_output[[epo]]=output
}

# Calculate the power and standard error
## utilize four test statistic functions (PRC, ERC, F_Sigma, GLR-l1) and take M1P1 and Bonf as baseline methods.

summ=read_experiments_epo(param_Hub_template,all_output)
#Power
power.table=extract_metric(metric_name='power')[c(1:8,14,17,22)]
colnames(power.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_sum","GLR_lasso")
#Standard error
SE.table=extract_metric(metric_name='SE')[c(1:8,14,17,22)]
colnames(SE.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_sum","GLR_lasso")



############### ER random graph (Q=c(0.2,0.4),Q_0=0.08,p=c(40,120),s=c(0.01,0.02),n=c(50,100)) ##########

param_ER_template <- list(
  model = "ER",
  N = c(50,100),
  P = c(40,120),
  modelParam = list(
    Q = c(0.2,0.4),
    Q0 = 0.08,
    s = c(0.01,0.02)
  )
)
M = 100
L = 3

path_prefix=paste0('Experiments-power-ER','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# repeat experiment 400 times and save the result
for (epo in 1:400){
  output = run_experiments_epo(param_ER_template,epo, M,L,List_CSS_stat =NULL ,dir_prefix=path_prefix)
  save(output,param_ER_template,file=paste0(path_prefix,epo,'.RData'))
}

#load the result
all_output=list()
for(epo in 1:400){
  load(file=paste0('Experiments-power-ER','/',epo,'.RData'))
  all_output[[epo]]=output
}

# Calculate the power and standard error
## utilize four test statistic functions (PRC, ERC, F_Sigma, GLR-l1) and take M1P1 and Bonf as baseline methods.

summ=read_experiments_epo(param_ER_template,all_output)
#Power
power.table=extract_metric(metric_name='power')[c(1:8,14,17,22)]
colnames(power.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_sum","GLR_lasso")
#Standard error
SE.table=extract_metric(metric_name='SE')[c(1:8,14,17,22)]
colnames(SE.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_sum","GLR_lasso")


########################################################################################################
    ###### Additional simulations on goodness-of-fit tests for Gaussian graphical models ######
########################################################################################################

## To test the GoF of GGMs, we additionally consider an experimental setup in Verzelen and Villers [2009]
## for a balanced comparison with existing methods

param_VV_template <- list(model = "Verzelen",
                          N = c(10,15,30),
                          P = c(15),
                          modelParam = list(
                            Q = c(0.1,0.4,1),
                            E = c(0.1,0.15)
                          )
)

path_prefix=paste0('Experiments-Verzelen','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# repeat experiment 400 times and save the result
for (epo in 1:400){
  output = run_experiments_epo(param_VV_template,epo, M,L,List_CSS_stat =NULL ,dir_prefix=path_prefix)
  save(output,param_VV_template,file=paste0(path_prefix,epo,'.RData'))
}

#load the result
all_output=list()
for(epo in 1:400){
  load(file=paste0('Experiments-Verzelen','/',epo,'.RData'))
  all_output[[epo]]=output
}

# Calculate the power and standard error
## utilize four test statistic functions (PRC, ERC, F_Sigma, F_max) and take M1P1 and Bonf as baseline methods.

summ=read_experiments_epo(param_ER_template,all_output)
#Power
power.table=extract_metric(metric_name='power')[c(1:8,14,17,16)]
colnames(power.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_sum","F_max")
#Standard error
SE.table=extract_metric(metric_name='SE')[c(1:8,14,17,16)]
colnames(SE.table) = c("n","p","K0","K","s","VV","DP","PRC","ERC","F_sum","F_max")


########################################################################################################
              ###### Additional simulations on GoF tests with prior information ######
########################################################################################################

# add weighting factor w=0.8
All_Stat_RC=function(x,graph,residual_type="PR",
                     threshold=0.1){
  p=ncol(x)
  if(residual_type=='PR'){
    pv_list=PRC_pvalues(x, graph)
  }else if(residual_type=='ER'){
    pv_list=ERC_pvalues(x,graph)
  }else {
    stop("Unknown residual type.")
  }
  
  all_rc_stat=c()
  
  for(i in 1:2){
    pv=pv_list[[i]] 
    
    z_score = - qnorm(pv/2)
    z_score = z_score*0.8
    for(moment_order in 2:1){
      # Calculate the final statistics based on the threshold
      s = 0
      for (i in 1:(p-1)){
        for (j in (i+1):p){
          if ( graph[i,j]==0 & pv[i, j] < threshold){
            s = s + abs(z_score[i, j])^moment_order
          }
        }
      }
      all_rc_stat=c(all_rc_stat,s)
    }
  }
  names(all_rc_stat)=paste0(ifelse(residual_type=='PR',"P","E"), 
                            c("RC_SS","RC_SA", "RC_Z_SS", "RC_Z_SA"))
  
  return(all_rc_stat)
}

############### band graph (K=6,K_0=1,p=c(20,120),s=0.1,n=c(20,40,80)) ##########

param_Band_template <- list(
  model = "Band_Band",
  N = c(20,40,80),
  P = c(20,120),
  modelParam = list(
    K0 = 1,
    K = 6,
    s = 0.1
  )
)
M = 100
L = 3

path_prefix=paste0('Experiments-prior','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# repeat experiment 400 times and save the result
for (epo in 1:400){
  output = run_experiments_epo(param_Band_template,epo, M,L,List_CSS_stat =NULL ,dir_prefix=path_prefix)
  save(output,param_Band_template,file=paste0(path_prefix,epo,'.RData'))
}

#load the result
all_output=list()
for(epo in 1:400){
  load(file=paste0('Experiments-prior','/',epo,'.RData'))
  all_output[[epo]]=output
}

# Calculate the power and standard error
## utilize four test statistic functions PRC-w and ERC-w.

summ=read_experiments_epo(param_Band_template,all_output)
#Power
power.table=extract_metric(metric_name='power')[c(1:5,8,14)]
colnames(power.table) = c("n","p","K0","K","s","PRC-w","ERC-w")
#Standard error
SE.table=extract_metric(metric_name='SE')[c(1:5,8,14)]
colnames(SE.table) = c("n","p","K0","K","s","PRC-w","ERC-w")
