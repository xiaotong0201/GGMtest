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
                                          # prop.test( mean(x < sig.lv), length(x), p=0.05, alternative='less',
                                          #           correct=F)$p.value
                                          # if(mean(vv<sig.lv)>0){
                                          # t.test(  c(vv<sig.lv) - sig.lv, alternative='greater' )$p.value
                                          # }else{NA}
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



if(1==2){


ind.experiment=0
all_output=list()
for(epo in 1:1200){
  load(file=paste0('Experiment_Results/Experiments-',ind.experiment,'/',epo,'.RData'))
  # load(file=paste0('../../../Programming/LocalExperiments/GGM/GoF/Experiments-',ind.experiment,'/',epo,'.RData'))
  all_output[[epo]]=output
}

summ=read_experiments_epo(param_template,all_output)


power.table=extract_metric('power')
SE.table=extract_metric('SE')
ospt.table=extract_metric('ospt')

sig.lv=0.05



one.sample.prop.test


min(one.sample.prop.test, na.rm=T)*6

min(      one.sample.prop.test[c("PRC_SS","ERC_Z_SS", "F_max", "F_sum","glr_glasso"),]  
    , na.rm=T)*6


options(digits=3)

# sort.table=power.table %>% arrange(E,Q) # verzelen
# subset(sort.table[,c(1:6,7,13,15,16)], E==0.1)
# subset(sort.table[,c(1:6,7,13,15,16)], E==0.15)


sort.table=power.table %>% arrange(p)
  print(sort.table)
print(SE.table %>% arrange(p))


### for analysis of null
if(ind.experiment==0){
  
ratio= (power.table-0.05)/SE.table


ratio=unlist((power.table[,c("PRC_SS","ERC_Z_SS", "F_max", "F_sum","glr_glasso")]-0.05)/
  SE.table[,c("PRC_SS","ERC_Z_SS", "F_max", "F_sum","glr_glasso")])
ratio[is.infinite(ratio)]=0

which.max(ratio)

hist(ratio,breaks=20)



pvs=summ[[3]]$pvalue.table[,'PRC_SS']

pvs=summ[[5]]$pvalue.table[,'glr_glasso']
mean(pvs<0.1)

hist(pvs)
qqplot((1:epo-0.5)/epo, pvs,pch=18,  cex=0.2 )
qqline(pvs,  dist=qunif)
pnorm( max(ratio) , lower=F) * length(ratio)

pnorm( max(ratio) , lower=F) * length(summ)

library(ggplot2)

## need sqrt(0.05*0.95/1600) ~ 0.00545 to make sure accurate

# Create the QQ plot using ggplot2's stat_qq function
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


ggsave("qqplot-null-graph-glr.pdf", plot = qq_plot, device = "pdf", width = 6, height = 6)


}

Is.Abs.Win=function(v,gap=0.1){
  second.win=sort(v,decreasing = T)[2]
  return(v>second.win+gap)
}

Is.Near.Win=function(v, tol= 0.05 ){
  return( tol + v>= max(v) )
}

abs.winner.table=t(apply(power.table[,-(1:5)],1, Is.Abs.Win))
colSums(abs.winner.table)



near.winner.table=t(apply(power.table[,-(1:5)],1,Is.Near.Win ))
colSums(near.winner.table)



power.table[which(power.table$VV > power.table$PRC_SS),]   # VV may works only in low dimensional (band, Hub, not ER)
power.table[which(power.table$DP > power.table$PRC_SS),]  # similarly



CONST_TOLERANCE=1/17

## compare SA and SS
## we can safely focus on SS only
power.table[which(power.table$PRC_SA > power.table$PRC_SS  + CONST_TOLERANCE ),]  # 
power.table[which(power.table$PRC_SA < power.table$PRC_SS - CONST_TOLERANCE),]  # 
power.table[which(power.table$ERC_SA > power.table$ERC_SS  + CONST_TOLERANCE ),]  # 
power.table[which(power.table$ERC_SA < power.table$ERC_SS - CONST_TOLERANCE),]  # ERCSA is much lower power 

power.table[which(power.table$PRC_Z_SA > power.table$PRC_Z_SS  + CONST_TOLERANCE ),]  # 
power.table[which(power.table$PRC_Z_SA < power.table$PRC_Z_SS - CONST_TOLERANCE),]  # 
power.table[which(power.table$ERC_Z_SA > power.table$ERC_Z_SS  + CONST_TOLERANCE ),]  # 
power.table[which(power.table$ERC_Z_SA < power.table$ERC_Z_SS - CONST_TOLERANCE),]  # ERCSA is much lower power 



## compare Z-stablizing
power.table[which(power.table$PRC_Z_SS > power.table$PRC_SS + CONST_TOLERANCE ),]  # 
power.table[which(power.table$PRC_Z_SS < power.table$PRC_SS - CONST_TOLERANCE ),]  # 
power.table[which(power.table$PRC_Z_SA > power.table$PRC_SA + CONST_TOLERANCE ),]  # 
power.table[which(power.table$PRC_Z_SA < power.table$PRC_SA - CONST_TOLERANCE ),]  # 

power.table[which(power.table$ERC_Z_SS > power.table$ERC_SS + CONST_TOLERANCE ),]  # E-Z better in lower dimension
power.table[which(power.table$ERC_Z_SS < power.table$ERC_SS - CONST_TOLERANCE),]  # E-Z worse in high-d
power.table[which(power.table$ERC_Z_SA > power.table$ERC_SA + CONST_TOLERANCE ),]  # 
power.table[which(power.table$ERC_Z_SA < power.table$ERC_SA - CONST_TOLERANCE),]  # E-Z worse in high-d


## compare F-stat
power.table[which(power.table$F_Z > power.table$PRC_SS + CONST_TOLERANCE ),] # F_Z is more powerful in high-dimension and sparse
power.table[which(power.table$F_Z < power.table$PRC_SS - CONST_TOLERANCE),] # 


power.table[which(power.table$F_sum > power.table$PRC_SS + CONST_TOLERANCE ),] # F_sum also works more powerful in high-dimension
power.table[which(power.table$F_sum < power.table$PRC_SS - CONST_TOLERANCE),] #  
power.table[which(power.table$F_sum > power.table$F_Z + CONST_TOLERANCE ),] # F_sum dominated by F_Z
power.table[which(power.table$F_sum < power.table$F_Z - CONST_TOLERANCE),] #  


power.table[which(power.table$glr_ridge> power.table$PRC_SS + CONST_TOLERANCE ),]  # when low dimension, glr_ridge works
power.table[which(power.table$glr_ridge < power.table$PRC_SS -  CONST_TOLERANCE ),] #  glr_ridge fails in high-dimension


power.table[which(power.table$glr_glasso> power.table$PRC_SS  + CONST_TOLERANCE ),]  # glr_lasso works when true is strong and sparse ()
power.table[which(power.table$glr_glasso< power.table$PRC_SS  - CONST_TOLERANCE ),]  # glr_lasso fails when true is weak and diffuse 

power.table[which(power.table$glr_glasso> power.table$F_Z  + CONST_TOLERANCE ),]  # glr_lasso dominated by F_Z
power.table[which(power.table$glr_glasso< power.table$F_Z  - CONST_TOLERANCE ),]  # 

## Factor analysis
reduce=power.table[,1:5]
reduce$quantity=log((power.table$VV+power.table$DP)/power.table$PRC_SS)
anova.model=lm(quantity~.+.^2,data = reduce)
anova(anova.model)

## Chosen cases
### band graph
df=arrange(subset(x = power.table, p==120&s==0.15), K)
subset(x = power.table, p==120&s==0.15&K==2) # F_sum, F_Z best
subset(x = power.table, p==120&s==0.25&K==2) # F_sum, F_Z best
subset(x = power.table, p==120&s==0.15&K==6) # PRC best
subset(x = power.table, p==120&s==0.2&K==6) # every one powerful, except ERC_SS ERC_SA

## Hub graph
round(arrange(subset(x = power.table[,c(1:7,8,12,16,18,19,20)], p==120), noise),4) #glr glasso especially well
round(arrange(subset(x = power.table[,c(1:7,8,12,16,18,19,20)], p==20), noise),4)

## ER graph
arrange(subset(x = power.table, p==40&Q==0.25), s)
arrange(subset(x = power.table, p==120&Q==0.25), s)

# arrange(subset(x = power.table, p==40&Q==0.5), s)
# arrange(subset(x = power.table, p==120&Q==0.5), s)


########################
}