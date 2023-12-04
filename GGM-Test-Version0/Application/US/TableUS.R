library(dplyr)
library(tables)
source("~/Library/CloudStorage/GoogleDrive-statdongminghuang@gmail.com/My Drive/Codes in R/ArrayToTable-Package.R")

################## Create Table

load('ExperimentResults/US-subsampling.RData')

sig.lv=0.05
power_table=c()
SE_table=c()
for(i in 1:length(subsampling_results)){
  Rejection=t(sapply(subsampling_results[[i]],function(x){c(x[1:2], unlist(x[-(1:2)]))})) < sig.lv
  power_table=rbind(power_table,
                    colMeans(Rejection)
  )
  
  SE_table=rbind(SE_table,
                 apply(Rejection,
                       2,function(v)sqrt(var(v)/length(v)))
  )
}

colnames(power_table) = colnames(SE_table)=c("VV","DP","PRC_SS","ERC_Z_SS", 'F_max',"F_sum","glr_glasso")

power_table=power_table[-5,]
SE_table=SE_table[-5,]

power_table=cbind(subsample_size=sample_sizes, power_table)
SE_table=cbind(subsample_size=sample_sizes, SE_table)


df_power=data.frame(power_table)
df_power_long <- tidyr::pivot_longer(df_power, cols = - 1  , 
                                     names_to = "Method", values_to = "Power")

df_SE=data.frame(SE_table)
df_SE_long <- tidyr::pivot_longer(df_SE, cols = - 1  ,
                                  names_to = "Method", values_to = "SE")

df_merged <- merge(df_power_long, df_SE_long, by = c( 'subsample_size', "Method"))

df_merged$subsample_size = as.factor(df_merged$subsample_size)
df_merged$Method = as.factor(df_merged$Method)


selected_methods=c("VV","DP","PRC_SS","ERC_Z_SS","F_sum","glr_glasso")

df_display=df_merged
df_display=subset(df_display,Method%in%selected_methods)
df_display$Method=factor(df_display$Method,levels=selected_methods)


df_display <- df_display %>% 
  group_by(across(1)) %>% 
  mutate(Is_Max = ifelse(Power == max(Power), TRUE, FALSE)) %>%
  ungroup()

df_display$Power<- sprintf("%.3f",df_display$Power)
df_display$Power<- ifelse(df_display$Is_Max, 
                          paste0("\\textbf{", df_display$Power, "}"),
                          df_display$Power)

df_display$Display <- sprintf("%s (%s)", 
                              df_display$Power, 
                              sub("^0", "", sprintf("%.3f",df_display$SE))
)              


tab <- tabular( subsample_size  ~ Method * Display  *c , data = df_display)
rlt=rowLabels(tab)
colnames(rlt)=c('$n_{s}$')
rowLabels(tab)<-rlt


clt=colLabels(tab)
clt  = clt[1:2,]
clt[2,]=c('VV','DP','PRC','ERC', 'F$_{\\Sigma}$',"GLR-${\\ell_1}$" )
colLabels(tab)<-clt

latex.tabular(tab,file ='US_subsample_power.tex', sepR=4, sepC=6)
