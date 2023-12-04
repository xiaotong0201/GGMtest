library(dplyr)

## For anaylsis
read_experiments_epo <- function(param_template, output_list,sig.lv=0.05) {
  MaxEpo=length(output_list)
  
  summary_result <- list()
  lFR=0
  
  for (ind.N in 1:length(param_template$N)) {
    N=param_template$N[ind.N]
    for (ind.P in 1:length(param_template$P)) {
      P=param_template$P[ind.P]
      
      combos <- expand.grid(param_template$modelParam)
      
      # Loop over each combination of model parameters
      for (ind.setup in 1:nrow(combos)) {
        
        
        np.id=(ind.N-1)*length(param_template$P)+ind.P
        exp.id=ind.setup+nrow(combos)*(np.id-1)
        
        pvalue.table=c()
        for(epo in 1:MaxEpo){
          pvalue.table= rbind(pvalue.table, c(unlist(output_list[[epo]][[exp.id]]) ))
        }
        
        colnames(pvalue.table)=c("VV","DP",names(output_list[[1]][[1]]$CSS_test_pvalues))
        
        # colMeans(pvalue.table)
        power=colMeans(pvalue.table < sig.lv)
        SE=apply(pvalue.table < sig.lv,2,function(v)sqrt(var(v)/length(v)))
        
        # print(power)
        
        modelParam = as.list(combos[ind.setup, ])
        names(modelParam)=names(param_template$modelParam)

        summary_result[[lFR+1]]  <- 
          list(model=param_template$model, 
               N=N,P = P, modelParam=modelParam,
               pvalue.table=pvalue.table,
               power=power,SE=SE
               )
        lFR=lFR+1
        
      }
    }
  }
  return(summary_result)
}


extract_metric=function(metric_name, sig.lv=0.05){
  if(metric_name=='ospt'){ # one sample proportional test
    one.sample.prop.test=sapply(summ, 
                                function(res){
                                  apply(res$pvalue.table,2, 
                                        function(vv){
                                          prop.test( sum(vv < sig.lv), length(vv), p=0.05, alternative='greater',correct=F)$p.value
                                          
                                        }
                                  )
                                })
    new_table=sapply(summ,function(x)(c(x$N,x$P,unlist(x$modelParam),apply(x$pvalue.table,2, 
                                                                           function(vv){
                                                                             prop.test( sum(vv < sig.lv), length(vv), p=0.05, alternative='greater',correct=F)$p.value
                                                                           })
    )))
  }else{
    new_table=sapply(summ,function(x)(c(x$N,x$P,unlist(x$modelParam),x[[metric_name]])))
  }
  rownames(new_table)[1:2]=c('n','p')
  new_table=t(new_table)
  new_table=as.data.frame(new_table)
  return(new_table)
}

