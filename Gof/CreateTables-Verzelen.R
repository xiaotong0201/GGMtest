source('PostAnalysis.R')
library(tables)
source("~/Library/CloudStorage/GoogleDrive-statdongminghuang@gmail.com/My Drive/Codes in R/ArrayToTable-Package.R")

ind.experiment=4

ETA=c(0.1,0.15)
ind.eta=1

for(ind.eta in 1:2){
    
    all_output=list()
    for(epo in 1:400){
      load(file=paste0('../../../Programming/LocalExperiments/GGM/GoF/Experiments-',ind.experiment,'/',epo,'.RData'))
      all_output[[epo]]=output
    }
    
    summ=read_experiments_epo(param_template,all_output)
    
    extract_metric=function(metric_name){
      new_table=sapply(summ,function(x)(c(x$N,x$P,unlist(x$modelParam),x[[metric_name]])))
      rownames(new_table)[1:2]=c('n','p')
      new_table=t(new_table)
      new_table=as.data.frame(new_table)
      return(new_table)
    }
    
    power.table=extract_metric('power')
    SE.table=extract_metric('SE')
    
    const.E=ETA[ ind.eta]
    
      
      df_power=subset(x = power.table, E==const.E)
      df_SE=subset(x = SE.table, E==const.E)

    
    num_factors=4
    
    # Assuming your data frame is named 'df'
    df_power_long <- tidyr::pivot_longer(df_power, cols = - (1:num_factors)  , 
                                         names_to = "Method", values_to = "Power")
    
    df_SE_long <- tidyr::pivot_longer(df_SE, cols = - (1:num_factors)  ,
                                      names_to = "Method", values_to = "SE")
    
    df_merged <- merge(df_power_long, df_SE_long, by = c( colnames(df_power)[1:num_factors], "Method"))
    for(j in 1:(num_factors+1)){
      df_merged[[j]]=as.factor(df_merged[[j]])
    }
    
    
    
    selected_methods=c("VV","DP","PRC_SS","ERC_Z_SS","F_sum","F_max")
    
    df_display=subset(df_merged,Method%in%selected_methods)
    df_display$Method=factor(df_display$Method,levels=selected_methods)
    
    
    
    
    
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
    


    
      
      tab <- tabular(Q  * n  ~ Method * Display  *c , data = df_display)
      rlt=rowLabels(tab)
      colnames(rlt)=c('$q$','$n$')
      
      
    rowLabels(tab)<-rlt
    
    
    clt=colLabels(tab)
    clt  = clt[1:2,]
    
    clt[2,]=c('VV','DP','PRC','ERC', 'F$_{\\Sigma}$',"F-max" )
  
    colLabels(tab)<-clt
    
    
    latex.tabular(tab, file = paste0('table-',ind.experiment,'-eta=',ETA[ind.eta],'.tex'),
                  sepR =3,
                  sepC=6)
    

}
