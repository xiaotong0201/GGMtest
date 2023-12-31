library(dplyr)
source('PostAnalysis.R')
library(tables)
source("~/Library/CloudStorage/GoogleDrive-statdongminghuang@gmail.com/My Drive/Codes in R/ArrayToTable-Package.R")

ind.experiment=3

const.K=6
# const.K=2

for(ind.experiment in 0:3){
  
  Max_epo=ifelse(ind.experiment==0,1200, 400)
  
  
  all_output=list()
  for(epo in 1:Max_epo){
    load(file=paste0('Experiment_Results/Experiments-',ind.experiment,'/',epo,'.RData'))
    all_output[[epo]]=output
  }
  
  summ=read_experiments_epo(param_template,all_output)
  
  power.table=extract_metric('power')
  SE.table=extract_metric('SE')
  
  
    if(ind.experiment==1){
      df_power=subset(x = power.table, K==const.K)
      df_SE=subset(x = SE.table,  K==const.K)
    }else if(ind.experiment==2){
      
      
      df_power=subset(x = power.table,  noise!=0.5)
      df_SE=subset(x = SE.table,  noise!=0.5)
      
    }else if(ind.experiment==3){
      
      df_power=power.table
      df_SE=SE.table 
    }else if(ind.experiment==0){
      
      df_power=power.table
      df_SE = extract_metric('ospt') # put the p-value in parthesis
    }
    
    
    num_factors=5
    
    # Assuming your data frame is named 'df'
    df_power_long <- tidyr::pivot_longer(df_power, cols = - (1:num_factors)  , 
                                         names_to = "Method", values_to = "Power")
    
    df_SE_long <- tidyr::pivot_longer(df_SE, cols = - (1:num_factors)  ,
                                      names_to = "Method", values_to = "SE")
    
    df_merged <- merge(df_power_long, df_SE_long, by = c( colnames(power.table)[1:num_factors], "Method"))
    
    
    
    for(j in 1:(num_factors+1)){
      df_merged[[j]]=as.factor(df_merged[[j]])
    }
    
    if(ind.experiment==2){
      df_merged$noise=factor(df_merged$noise, levels=c(2,1.4,0.9))
    }
    
    
    selected_methods=c("VV","DP","PRC_SS","ERC_Z_SS","F_sum","glr_glasso")
    
    df_display=df_merged
    df_display=subset(df_display,Method%in%selected_methods)
    df_display$Method=factor(df_display$Method,levels=selected_methods)
    
    
    # ## highlight the combined 
    # df_display$Display <- sprintf("%.3f (%s)", 
    #                               df_display$Power, 
    #                               sub("^0", "", sprintf("%.3f",df_display$SE))
    #                               )
    # df_display <- df_display %>% 
    #   group_by(across(1:num_factors)) %>% 
    #   mutate(Is_Max = ifelse(Power == max(Power), TRUE, FALSE)) %>%
    #   ungroup()
    # ## highlight the combined 
    # df_display$Display <- ifelse(df_display$Is_Max, 
    #                              paste0("\\textbf{", df_display$Display, "}"), 
    #                              df_display$Display)
    
    
    if(ind.experiment!=0){
      df_display <- df_display %>% 
        group_by(across(1:num_factors)) %>% 
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
    }else{
      df_display$Display <- sprintf("%.3f (%s)",
                                    df_display$Power,
                                    sprintf("%.3f",df_display$SE)
      )
    }
    
    
    
    
    # ## highlight the combined 
    # df_display$Display <- ifelse(df_display$Is_Max, 
    #                              paste0("\\textbf{", df_display$Display, "}"), 
    #                              df_display$Display)
    
    
    
    if(ind.experiment==1){
      tab <- tabular(p* s * n  ~ Method * Display  *c , data = df_display)
      rlt=rowLabels(tab)
      colnames(rlt)=c('$p$','$s$','$n$')
    }else if(ind.experiment==2){
      tab <- tabular(p* noise * n  ~ Method * Display  *c , data = df_display)
      rlt=rowLabels(tab)
      colnames(rlt)=c('$p$','$\\xi$','$n$')
    }else if(ind.experiment==3){
      
      tab <- tabular(p* Q * s * n  ~ Method * Display  *c , data = df_display)
      rlt=rowLabels(tab)
      colnames(rlt)=c('$p$', '$q$','$s$','$n$')
    }else if(ind.experiment==0){
      
      tab <- tabular(p*s * n  ~ Method * Display  *c , data = df_display)
      rlt=rowLabels(tab)
      colnames(rlt)=c('$p$',  '$s$','$n$')
      
    }
    
    rowLabels(tab)<-rlt
    
    
    clt=colLabels(tab)
    clt  = clt[1:2,]
    # clt[2,]=c('VV','DP','PRC$_{\\mbox{SS}}$','ERC$_{Z}$', 'F$_{\\Sigma}$',"GLR$_{\\ell 1}$" )
    
    clt[2,]=c('VV','DP','PRC','ERC', 'F$_{\\Sigma}$',"GLR-${\\ell_1}$" )
    colLabels(tab)<-clt
    
    
    latex.tabular(tab, file = paste0('table-',ind.experiment,
                                     ifelse(ind.experiment==1, paste0('-K=',const.K), ''), 
                                     '.tex'),
                  sepR = c(3, 3,3,2)[ind.experiment+1]  ,
                  sepC=6   )
    
  
  
  
  
}
