
### We summariz the setup for generating data and running CRT test of the simulation studies

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
                ###### Gaussian linear regression ######
##########################################################################

### Low-dimensional

# Set p = 20, N = 50. 
# We examine our G-CRT procedure with statistics LM-SST and LM-SSR by comparison  to the 
# classical F-test and the dCRT using a Gaussian Lasso model (Bonferroni adjusted)

# Set parameters
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

# Set the storage path
path_prefix=paste0('Experiments-llw','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# Repeat experiment and save results
for (epo in 1:400){
  output = CRT(popul_param_llw,epo = epo)
  save(output,popul_param_llw,file=paste0(path_prefix,epo,'.RData'))
}

# Load results
all_output=list()
for(epo in 1:400){
  load(file=paste0(path_prefix,epo,'.RData'))
  all_output[[epo]]=output
}

combined_pvalues=list()
for(ind.signal in 1:length(popul_param_llw$modelParam$grid_signal)){
  combined_pvalues[[ind.signal]]=t(sapply(all_output, function(v)v[[ind.signal]]$pvalues))
}
combined_pvalues=lapply(combined_pvalues,function(x)x[,-c(5)])
combined_pvalues=lapply(combined_pvalues,function(x){
  colnames(x)=c("LM-SST","LM-SSR","F-test","dCRT")
  return(x)})

#Present the power curve of various testing methods in the experiment, in which power is plotted against varying values of theta

Fun_SE=function(x){sd(x)/sqrt(length(x))}

pvalues_array=simplify2array(combined_pvalues)
sig.lv=0.05

power_table=data.frame(apply(pvalues_array<sig.lv, 3:2,mean)); 
power_SE_table=data.frame(apply(pvalues_array<sig.lv, 3:2,Fun_SE))

num_methods=ncol(power_table)
power_table$theta=popul_param_llw$modelParam$grid_signal
power_SE_table$theta=popul_param_llw$modelParam$grid_signal



df_melt=melt(power_table,id.vars="theta", variable_name = "Test")
df_melt$Power=df_melt$value;df_melt$value=NULL

SE_melt <- melt(power_SE_table, id.vars = "theta", variable_name = "Test")
SE_melt$SE=SE_melt$value;SE_melt$value=NULL

# Merge the data and SE data frames
merged_df <- merge(df_melt, SE_melt, by = c("theta", "Test"))

Select_Methods=c()
Select_Theta=c()

distinct_colors <- scales::hue_pal()(num_methods)
distinct_linetypes <- rep(c("solid", "dashed", "dotted", "twodash", "longdash", "dotdash"), length.out = num_methods)
distinct_shapes <- rep( c(1,2,3,4,5,6,8,15,16,17,18), length.out= num_methods)


ggplot(merged_df, aes(x = theta, y = Power, color = Test, shape=Test, linetype = Test)) +
  geom_point() +
  geom_line(size = 0.6) +
  geom_errorbar(aes(ymin=Power-SE, ymax=Power+SE), width=0.1) +
  # geom_ribbon(aes(ymin = Power - SE, ymax = Power + SE, fill = Test), alpha = 0.2) +
  theme(legend.position = "right") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype="dotted", color = "black") +
  annotate("text", x = max(merged_df$theta), y = 0.05, label = "Level", vjust = -1) +
  labs( #title = "Power Curve",
    x = expression(theta),
    y = "Power")  +
  scale_color_manual(values = distinct_colors) +
  scale_shape_manual(values = distinct_shapes)+
  scale_fill_manual(values = distinct_colors) 

### High-dimensional

# Set p = 120, N = 80. 
# We examine our G-CRT procedure with statistics L1-R-SST and L1-R-SSR by comparison to 
# the Bonferroni adjusted tests using dCRT and the de-sparsified Lasso labeled as De-Lasso
 
# Set parameters
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

# Set the storage path
path_prefix=paste0('Experiments-hlw','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# Repeat experiment and save results
for (epo in 1:400){
  output = CRT(popul_param_hlw,epo = epo)
  save(output,popul_param_hlw,file=paste0(path_prefix,epo,'.RData'))
}

# Load results
all_output=list()
for(epo in 1:400){
  load(file=paste0(path_prefix,epo,'.RData'))
  all_output[[epo]]=output
}

combined_pvalues=list()
for(ind.signal in 1:length(popul_param_hlw$modelParam$grid_signal)){
  combined_pvalues[[ind.signal]]=t(sapply(all_output, function(v)v[[ind.signal]]$pvalues))
}
combined_pvalues=lapply(combined_pvalues,function(x)x[,-c(4,6)])
combined_pvalues=lapply(combined_pvalues,function(x){
  colnames(x)=c("LM-L1-R-SST","LM-L1-R-SSR","De-Lasso","dCRT")
  return(x)})

# Present the power curve of various testing methods in the experiment, in which power is plotted against varying values of theta

Fun_SE=function(x){sd(x)/sqrt(length(x))}

pvalues_array=simplify2array(combined_pvalues)
sig.lv=0.05

power_table=data.frame(apply(pvalues_array<sig.lv, 3:2,mean)); 
power_SE_table=data.frame(apply(pvalues_array<sig.lv, 3:2,Fun_SE))

num_methods=ncol(power_table)
power_table$theta=popul_param_hlw$modelParam$grid_signal
power_SE_table$theta=popul_param_hlw$modelParam$grid_signal



df_melt=melt(power_table,id.vars="theta", variable_name = "Test")
df_melt$Power=df_melt$value;df_melt$value=NULL

SE_melt <- melt(power_SE_table, id.vars = "theta", variable_name = "Test")
SE_melt$SE=SE_melt$value;SE_melt$value=NULL

# Merge the data and SE data frames
merged_df <- merge(df_melt, SE_melt, by = c("theta", "Test"))

Select_Methods=c()
Select_Theta=c()

distinct_colors <- scales::hue_pal()(num_methods)
distinct_linetypes <- rep(c("solid", "dashed", "dotted", "twodash", "longdash", "dotdash"), length.out = num_methods)
distinct_shapes <- rep( c(1,2,3,4,5,6,8,15,16,17,18), length.out= num_methods)


ggplot(merged_df, aes(x = theta, y = Power, color = Test, shape=Test, linetype = Test)) +
  geom_point() +
  geom_line(size = 0.6) +
  geom_errorbar(aes(ymin=Power-SE, ymax=Power+SE), width=0.1) +
  # geom_ribbon(aes(ymin = Power - SE, ymax = Power + SE, fill = Test), alpha = 0.2) +
  theme(legend.position = "right") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype="dotted", color = "black") +
  annotate("text", x = max(merged_df$theta), y = 0.05, label = "Level", vjust = -1) +
  labs( #title = "Power Curve",
    x = expression(theta),
    y = "Power")  +
  scale_color_manual(values = distinct_colors) +
  scale_shape_manual(values = distinct_shapes)+
  scale_fill_manual(values = distinct_colors) 

##########################################################################
              ###### Logistic regression ######
##########################################################################

### Low-dimensional

# Set p = 20, N = 50. 
# we evaluate the power of our G-CRT with the GLM-Dev statistic and the sum of importance scores from random forests (RF) using the mean decrease in accuracy
# also consider the Bonferroni adjusted dCRT with the logistic Lasso model and the classical Chi-squared test

# Set parameters
popul_param_lgw <- list(
  setting = "l_gl_w",
  N = c(50),
  P = c(20),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal=seq(0,3,by=0.5)
  )
)
# Set the storage path
path_prefix=paste0('Experiments-lgw','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# Repeat experiment and save results
for (epo in 1:400){
  output = CRT(popul_param_lgw,epo = epo)
  save(output,popul_param_lgw,file=paste0(path_prefix,epo,'.RData'))
}

# Load results
all_output=list()
for(epo in 1:400){
  load(file=paste0(path_prefix,epo,'.RData'))
  all_output[[epo]]=output
}

combined_pvalues=list()
for(ind.signal in 1:length(popul_param_lgw$modelParam$grid_signal)){
  combined_pvalues[[ind.signal]]=t(sapply(all_output, function(v)v[[ind.signal]]$pvalues))
}
combined_pvalues=lapply(combined_pvalues,function(x)x[,-c(5)])
combined_pvalues=lapply(combined_pvalues,function(x){
  colnames(x)=c("GLM-Dev","RF","Chi-square","dCRT")
  return(x)})

# Present the power curve of various testing methods in the experiment, in which power is plotted against varying values of theta
Fun_SE=function(x){sd(x)/sqrt(length(x))}

pvalues_array=simplify2array(combined_pvalues)
sig.lv=0.05

power_table=data.frame(apply(pvalues_array<sig.lv, 3:2,mean)); 
power_SE_table=data.frame(apply(pvalues_array<sig.lv, 3:2,Fun_SE))

num_methods=ncol(power_table)
power_table$theta=popul_param_lgw$modelParam$grid_signal
power_SE_table$theta=popul_param_lgw$modelParam$grid_signal



df_melt=melt(power_table,id.vars="theta", variable_name = "Test")
df_melt$Power=df_melt$value;df_melt$value=NULL

SE_melt <- melt(power_SE_table, id.vars = "theta", variable_name = "Test")
SE_melt$SE=SE_melt$value;SE_melt$value=NULL

# Merge the data and SE data frames
merged_df <- merge(df_melt, SE_melt, by = c("theta", "Test"))

Select_Methods=c()
Select_Theta=c()

distinct_colors <- scales::hue_pal()(num_methods)
distinct_linetypes <- rep(c("solid", "dashed", "dotted", "twodash", "longdash", "dotdash"), length.out = num_methods)
distinct_shapes <- rep( c(1,2,3,4,5,6,8,15,16,17,18), length.out= num_methods)


ggplot(merged_df, aes(x = theta, y = Power, color = Test, shape=Test, linetype = Test)) +
  geom_point() +
  geom_line(size = 0.6) +
  geom_errorbar(aes(ymin=Power-SE, ymax=Power+SE), width=0.1) +
  # geom_ribbon(aes(ymin = Power - SE, ymax = Power + SE, fill = Test), alpha = 0.2) +
  theme(legend.position = "right") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype="dotted", color = "black") +
  annotate("text", x = max(merged_df$theta), y = 0.05, label = "Level", vjust = -1) +
  labs( #title = "Power Curve",
    x = expression(theta),
    y = "Power")  +
  scale_color_manual(values = distinct_colors) +
  scale_shape_manual(values = distinct_shapes)+
  scale_fill_manual(values = distinct_colors) 

### High-dimensional

# Set p = 120, N = 80
# we evaluate the power of our G-CRT with GLM−L1−D, GLM−L1−R−SST, RF, CGM and dCRT

# Set parameters
popul_param_hgw <- list(
  setting = "h_gl_w",
  N = c(80),
  P = c(120),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= 5  #seq(5,by=1)
  )
)
# Set the storage path
path_prefix=paste0('Experiments-hgw','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# Repeat experiment and save results
for (epo in 1:2){
  output = CRT(popul_param_hgw,epo = epo)
  save(output,popul_param_hgw,file=paste0(path_prefix,epo,'.RData'))
}

# Load results
all_output=list()
for(epo in 1:2){
  load(file=paste0(path_prefix,epo,'.RData'))
  all_output[[epo]]=output
}

combined_pvalues=list()
for(ind.signal in 1:length(popul_param_hgw$modelParam$grid_signal)){
  combined_pvalues[[ind.signal]]=t(sapply(all_output, function(v)v[[ind.signal]]$pvalues))
}
combined_pvalues=lapply(combined_pvalues,function(x)x[,-c(6)])
combined_pvalues=lapply(combined_pvalues,function(x){
  colnames(x)=c("GLM-L1-D","GLM-L1-R-SST","RF","CGM","dCRT")
  return(x)})

# Present the power curve of various testing methods in the experiment, in which power is plotted against varying values of theta
Fun_SE=function(x){sd(x)/sqrt(length(x))}

pvalues_array=simplify2array(combined_pvalues)
sig.lv=0.05

power_table=data.frame(apply(pvalues_array<sig.lv, 3:2,mean)); 
power_SE_table=data.frame(apply(pvalues_array<sig.lv, 3:2,Fun_SE))

num_methods=ncol(power_table)
power_table$theta=popul_param_hgw$modelParam$grid_signal
power_SE_table$theta=popul_param_hgw$modelParam$grid_signal



df_melt=melt(power_table,id.vars="theta", variable_name = "Test")
df_melt$Power=df_melt$value;df_melt$value=NULL

SE_melt <- melt(power_SE_table, id.vars = "theta", variable_name = "Test")
SE_melt$SE=SE_melt$value;SE_melt$value=NULL

# Merge the data and SE data frames
merged_df <- merge(df_melt, SE_melt, by = c("theta", "Test"))

Select_Methods=c()
Select_Theta=c()

distinct_colors <- scales::hue_pal()(num_methods)
distinct_linetypes <- rep(c("solid", "dashed", "dotted", "twodash", "longdash", "dotdash"), length.out = num_methods)
distinct_shapes <- rep( c(1,2,3,4,5,6,8,15,16,17,18), length.out= num_methods)


ggplot(merged_df, aes(x = theta, y = Power, color = Test, shape=Test, linetype = Test)) +
  geom_point() +
  geom_line(size = 0.6) +
  geom_errorbar(aes(ymin=Power-SE, ymax=Power+SE), width=0.1) +
  # geom_ribbon(aes(ymin = Power - SE, ymax = Power + SE, fill = Test), alpha = 0.2) +
  theme(legend.position = "right") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype="dotted", color = "black") +
  annotate("text", x = max(merged_df$theta), y = 0.05, label = "Level", vjust = -1) +
  labs( #title = "Power Curve",
    x = expression(theta),
    y = "Power")  +
  scale_color_manual(values = distinct_colors) +
  scale_shape_manual(values = distinct_shapes)+
  scale_fill_manual(values = distinct_colors) 

##########################################################################
              ###### Nonlinear regression ######
##########################################################################

### Low-dimensional

# Set p = 20, N = 50.
# we evaluate the power of our G-CRT with RF, RF-D, RF-RR, F-test, dCRT(RF)

# Set parameters
popul_param_llm <- list(
  setting = "l_l_m",
  N = c(50),
  P = c(20),
  setT = 1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal=seq(0,10,by=0.5)
  )
)

path_prefix=paste0('Experiments-llm','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# Repeat experiment and save results
for (epo in 1:2){
  output = CRT(popul_param_llm,epo = epo)
  save(output,popul_param_llm,file=paste0(path_prefix,epo,'.RData'))
}

# Load results
all_output=list()
for(epo in 1:2){
  load(file=paste0(path_prefix,epo,'.RData'))
  all_output[[epo]]=output
}

combined_pvalues=list()
for(ind.signal in 1:length(popul_param_llm$modelParam$grid_signal)){
  combined_pvalues[[ind.signal]]=t(sapply(all_output, function(v)v[[ind.signal]]$pvalues))
}
combined_pvalues=lapply(combined_pvalues,function(x)x[,-c(6)])
combined_pvalues=lapply(combined_pvalues,function(x){
  colnames(x)=c("RF","RF-D","RF-RR","MARS","F-test","dCRT(RF)")
  return(x)})

# Present the power curve of various testing methods in the experiment, in which power is plotted against varying values of theta

Fun_SE=function(x){sd(x)/sqrt(length(x))}

pvalues_array=simplify2array(combined_pvalues)
sig.lv=0.05

power_table=data.frame(apply(pvalues_array<sig.lv, 3:2,mean)); 
power_SE_table=data.frame(apply(pvalues_array<sig.lv, 3:2,Fun_SE))

num_methods=ncol(power_table)
power_table$theta=popul_param_llm$modelParam$grid_signal
power_SE_table$theta=popul_param_llm$modelParam$grid_signal



df_melt=melt(power_table,id.vars="theta", variable_name = "Test")
df_melt$Power=df_melt$value;df_melt$value=NULL

SE_melt <- melt(power_SE_table, id.vars = "theta", variable_name = "Test")
SE_melt$SE=SE_melt$value;SE_melt$value=NULL

# Merge the data and SE data frames
merged_df <- merge(df_melt, SE_melt, by = c("theta", "Test"))

Select_Methods=c()
Select_Theta=c()

distinct_colors <- scales::hue_pal()(num_methods)
distinct_linetypes <- rep(c("solid", "dashed", "dotted", "twodash", "longdash", "dotdash"), length.out = num_methods)
distinct_shapes <- rep( c(1,2,3,4,5,6,8,15,16,17,18), length.out= num_methods)


ggplot(merged_df, aes(x = theta, y = Power, color = Test, shape=Test, linetype = Test)) +
  geom_point() +
  geom_line(size = 0.6) +
  geom_errorbar(aes(ymin=Power-SE, ymax=Power+SE), width=0.1) +
  # geom_ribbon(aes(ymin = Power - SE, ymax = Power + SE, fill = Test), alpha = 0.2) +
  theme(legend.position = "right") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype="dotted", color = "black") +
  annotate("text", x = max(merged_df$theta), y = 0.05, label = "Level", vjust = -1) +
  labs( #title = "Power Curve",
    x = expression(theta),
    y = "Power")  +
  scale_color_manual(values = distinct_colors) +
  scale_shape_manual(values = distinct_shapes)+
  scale_fill_manual(values = distinct_colors) 

### High-dimensional

# Set p = 120, N = 80
# we evaluate the power of our G-CRT with RF, RF-D, RF-RR, De-Lasso, dCRT(RF)
# Set parameters
popul_param_hlm <- list(
  setting = "h_l_m",
  N = c(80),
  P = c(120),
  setT = 1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= 15 # seq(0,15,by=3)
  )
)

path_prefix=paste0('Experiments-hlm','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# Repeat experiment and save results
for (epo in 1:2){
  output = CRT(popul_param_hlm,epo = epo)
  save(output,popul_param_hlm,file=paste0(path_prefix,epo,'.RData'))
}

# Load results
all_output=list()
for(epo in 1:2){
  load(file=paste0(path_prefix,epo,'.RData'))
  all_output[[epo]]=output
}

combined_pvalues=list()
for(ind.signal in 1:length(popul_param_hlm$modelParam$grid_signal)){
  combined_pvalues[[ind.signal]]=t(sapply(all_output, function(v)v[[ind.signal]]$pvalues))
}
combined_pvalues=lapply(combined_pvalues,function(x)x[,-c(6,7)])
combined_pvalues=lapply(combined_pvalues,function(x){
  colnames(x)=c("RF","RF-D","RF-RR","MARS","De-Lasso","dCRT(RF)")
  return(x)})

# Present the power curve of various testing methods in the experiment, in which power is plotted against varying values of theta
Fun_SE=function(x){sd(x)/sqrt(length(x))}

pvalues_array=simplify2array(combined_pvalues)
sig.lv=0.05

power_table=data.frame(apply(pvalues_array<sig.lv, 3:2,mean)); 
power_SE_table=data.frame(apply(pvalues_array<sig.lv, 3:2,Fun_SE))

num_methods=ncol(power_table)
power_table$theta=popul_param_hlm$modelParam$grid_signal
power_SE_table$theta=popul_param_hlm$modelParam$grid_signal



df_melt=melt(power_table,id.vars="theta", variable_name = "Test")
df_melt$Power=df_melt$value;df_melt$value=NULL

SE_melt <- melt(power_SE_table, id.vars = "theta", variable_name = "Test")
SE_melt$SE=SE_melt$value;SE_melt$value=NULL

# Merge the data and SE data frames
merged_df <- merge(df_melt, SE_melt, by = c("theta", "Test"))

Select_Methods=c()
Select_Theta=c()

distinct_colors <- scales::hue_pal()(num_methods)
distinct_linetypes <- rep(c("solid", "dashed", "dotted", "twodash", "longdash", "dotdash"), length.out = num_methods)
distinct_shapes <- rep( c(1,2,3,4,5,6,8,15,16,17,18), length.out= num_methods)


ggplot(merged_df, aes(x = theta, y = Power, color = Test, shape=Test, linetype = Test)) +
  geom_point() +
  geom_line(size = 0.6) +
  geom_errorbar(aes(ymin=Power-SE, ymax=Power+SE), width=0.1) +
  # geom_ribbon(aes(ymin = Power - SE, ymax = Power + SE, fill = Test), alpha = 0.2) +
  theme(legend.position = "right") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype="dotted", color = "black") +
  annotate("text", x = max(merged_df$theta), y = 0.05, label = "Level", vjust = -1) +
  labs( #title = "Power Curve",
    x = expression(theta),
    y = "Power")  +
  scale_color_manual(values = distinct_colors) +
  scale_shape_manual(values = distinct_shapes)+
  scale_fill_manual(values = distinct_colors) 

##########################################################################
                ###### Nonlinear binary regression ######
##########################################################################

### Low-dimensional

# Set p = 20, N = 50
# we evaluate the power of our G-CRT with RF, RF-D, RF-RR, Chi−square, dCRT(RF)

# Set parameters
popul_param_lgm <- list(
  setting = "l_gl_m",
  N = c(50),
  P = c(20),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal=c(3,5)  #seq(0,3,by=0.5)
  )
)
#
path_prefix=paste0('Experiments-lgm','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# Repeat experiment and save results
for (epo in 1:2){
  output = CRT(popul_param_lgm,epo = epo)
  save(output,popul_param_lgm,file=paste0(path_prefix,epo,'.RData'))
}

# Load results
all_output=list()
for(epo in 1:2){
  load(file=paste0(path_prefix,epo,'.RData'))
  all_output[[epo]]=output
}

combined_pvalues=list()
for(ind.signal in 1:length(popul_param_lgm$modelParam$grid_signal)){
  combined_pvalues[[ind.signal]]=t(sapply(all_output, function(v)v[[ind.signal]]$pvalues))
}
combined_pvalues=lapply(combined_pvalues,function(x)x[,-c(5)])
combined_pvalues=lapply(combined_pvalues,function(x){
  colnames(x)=c("RF","RF-D","RF-RR","Chi-square","dCRT(RF)")
  return(x)})

# Present the power curve of various testing methods in the experiment, in which power is plotted against varying values of theta
Fun_SE=function(x){sd(x)/sqrt(length(x))}

pvalues_array=simplify2array(combined_pvalues)
sig.lv=0.05

power_table=data.frame(apply(pvalues_array<sig.lv, 3:2,mean)); 
power_SE_table=data.frame(apply(pvalues_array<sig.lv, 3:2,Fun_SE))

num_methods=ncol(power_table)
power_table$theta=popul_param_lgm$modelParam$grid_signal
power_SE_table$theta=popul_param_lgm$modelParam$grid_signal



df_melt=melt(power_table,id.vars="theta", variable_name = "Test")
df_melt$Power=df_melt$value;df_melt$value=NULL

SE_melt <- melt(power_SE_table, id.vars = "theta", variable_name = "Test")
SE_melt$SE=SE_melt$value;SE_melt$value=NULL

# Merge the data and SE data frames
merged_df <- merge(df_melt, SE_melt, by = c("theta", "Test"))

Select_Methods=c()
Select_Theta=c()

distinct_colors <- scales::hue_pal()(num_methods)
distinct_linetypes <- rep(c("solid", "dashed", "dotted", "twodash", "longdash", "dotdash"), length.out = num_methods)
distinct_shapes <- rep( c(1,2,3,4,5,6,8,15,16,17,18), length.out= num_methods)


ggplot(merged_df, aes(x = theta, y = Power, color = Test, shape=Test, linetype = Test)) +
  geom_point() +
  geom_line(size = 0.6) +
  geom_errorbar(aes(ymin=Power-SE, ymax=Power+SE), width=0.1) +
  # geom_ribbon(aes(ymin = Power - SE, ymax = Power + SE, fill = Test), alpha = 0.2) +
  theme(legend.position = "right") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype="dotted", color = "black") +
  annotate("text", x = max(merged_df$theta), y = 0.05, label = "Level", vjust = -1) +
  labs( #title = "Power Curve",
    x = expression(theta),
    y = "Power")  +
  scale_color_manual(values = distinct_colors) +
  scale_shape_manual(values = distinct_shapes)+
  scale_fill_manual(values = distinct_colors) 

### High-dimensional

# Set p = 120, N = 80. 
# we evaluate the power of our G-CRT with RF, RF-D, RF-RR, CGM, dCRT(RF)

# Set parameters
popul_param_hgm <- list(
  setting = "h_gl_m",
  N = c(80),
  P = c(120),
  setT=1:8, 
  modelParam = list(
    K = 6,
    s = 0.2,
    grid_signal= 5  # 0:5
  )
)
#
path_prefix=paste0('Experiments-hgm','/')
if(!dir.exists(path_prefix)){
  dir.create(path_prefix)
}

# Repeat experiment and save results
for (epo in 1:2){
  output = CRT(popul_param_hgm,epo = epo)
  save(output,popul_param_hgm,file=paste0(path_prefix,epo,'.RData'))
}

# Load results
all_output=list()
for(epo in 1:2){
  load(file=paste0(path_prefix,epo,'.RData'))
  all_output[[epo]]=output
}

combined_pvalues=list()
for(ind.signal in 1:length(popul_param_hgm$modelParam$grid_signal)){
  combined_pvalues[[ind.signal]]=t(sapply(all_output, function(v)v[[ind.signal]]$pvalues))
}
combined_pvalues=lapply(combined_pvalues,function(x)x[,-c(5)])
combined_pvalues=lapply(combined_pvalues,function(x){
  colnames(x)=c("RF","RF-D","RF-RR","CGM","dCRT(RF)")
  return(x)})

# Present the power curve of various testing methods in the experiment, in which power is plotted against varying values of theta
Fun_SE=function(x){sd(x)/sqrt(length(x))}

pvalues_array=simplify2array(combined_pvalues)
sig.lv=0.05

power_table=data.frame(apply(pvalues_array<sig.lv, 3:2,mean)); 
power_SE_table=data.frame(apply(pvalues_array<sig.lv, 3:2,Fun_SE))

num_methods=ncol(power_table)
power_table$theta=popul_param_hgm$modelParam$grid_signal
power_SE_table$theta=popul_param_hgm$modelParam$grid_signal



df_melt=melt(power_table,id.vars="theta", variable_name = "Test")
df_melt$Power=df_melt$value;df_melt$value=NULL

SE_melt <- melt(power_SE_table, id.vars = "theta", variable_name = "Test")
SE_melt$SE=SE_melt$value;SE_melt$value=NULL

# Merge the data and SE data frames
merged_df <- merge(df_melt, SE_melt, by = c("theta", "Test"))

Select_Methods=c()
Select_Theta=c()

distinct_colors <- scales::hue_pal()(num_methods)
distinct_linetypes <- rep(c("solid", "dashed", "dotted", "twodash", "longdash", "dotdash"), length.out = num_methods)
distinct_shapes <- rep( c(1,2,3,4,5,6,8,15,16,17,18), length.out= num_methods)


ggplot(merged_df, aes(x = theta, y = Power, color = Test, shape=Test, linetype = Test)) +
  geom_point() +
  geom_line(size = 0.6) +
  geom_errorbar(aes(ymin=Power-SE, ymax=Power+SE), width=0.1) +
  # geom_ribbon(aes(ymin = Power - SE, ymax = Power + SE, fill = Test), alpha = 0.2) +
  theme(legend.position = "right") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype="dotted", color = "black") +
  annotate("text", x = max(merged_df$theta), y = 0.05, label = "Level", vjust = -1) +
  labs( #title = "Power Curve",
    x = expression(theta),
    y = "Power")  +
  scale_color_manual(values = distinct_colors) +
  scale_shape_manual(values = distinct_shapes)+
  scale_fill_manual(values = distinct_colors) 
