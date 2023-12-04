library(dplyr)
library(tables)
source("~/Library/CloudStorage/GoogleDrive-statdongminghuang@gmail.com/My Drive/Codes in R/ArrayToTable-Package.R")

########### Graph G

G_pv_table=data.frame(
  Method=factor(c('VV','DP','Fsum'),levels=c('VV','DP','Fsum')),
  pValue=c(GoF_result[[1]], GoF_result[[2]], GoF_result[[3]][[1]])
)

G_pv_table=G_pv_table%>%
  mutate(pValue=sprintf("%.3f", pValue))


tab=tabular( pValue  ~ Method*c ,G_pv_table)

colLabels(tab)=colLabels(tab)[1:2,]
rowLabels(tab)[1,1]='P-Value'
colLabels(tab)[2,]=c('VV','DP','F$_{\\Sigma}$')
latex.tabular(tab,file='Stock_GoF.tex',sepC = 3, sepR = 1)

##############

########  Type 1 error ######



Fn_SE=function(v){sd(v)/sqrt(length(v))}

null_from_data_to_frame=function(null_result,char_setT){
  
  ind=c(1,2,3,4,5)
  
  
  null_pvalue_table = t(sapply(null_result,  '[[', "pvalues"))
  power=colMeans(null_pvalue_table < 0.05)
  se=apply(null_pvalue_table,2,Fn_SE)
temp_table=data.frame(
  Method=factor(c('LM-L1-R-SSR','RF-RR','De-Lasso', 'dCRT (L1)', 'dCRT (RF)'),
                levels = c('LM-L1-R-SSR','RF-RR','De-Lasso', 'dCRT (L1)', 'dCRT (RF)') ),
  power=power[ind], 
  se=se[ind]
)
temp_table=temp_table%>%
  mutate(setT=char_setT,Error=sprintf("%.3f (%.3f)", power, se))
return(temp_table)
}

load('Type-I_Weak5_DIA_2023-11-12.RData')
type1_table_weak=null_from_data_to_frame(null_result, char_setT="$\\mc{T}_1$")


load('Type-I_Strong4_DIA_2023-11-12.RData')
type1_table_strong=null_from_data_to_frame(null_result, char_setT="$\\mc{T}_2$")

type1error_table=rbind(type1_table_weak, type1_table_strong)

type1error_table$setT=as.factor(type1error_table$setT)

tab=tabular( setT ~ Method*Error*c ,type1error_table)


colLabels(tab)=colLabels(tab)[1:2,]
attr(rowLabels(tab),'dimnames')[[2]]='\\mc{T}'
latex.tabular(tab,file='Stock_Type-I-error.tex',sepC = 5, sepR = 2)





######### p value ##########






CI_from_data_to_frame=function(CI_test_result,    char_setT){
  
  ind=c(1,2,4,5)
  temp_table=data.frame(
    Method=factor(c('LM-L1-R-SSR','RF-RR', 'dCRT (L1)', 'dCRT (RF)'),
                  levels=c('LM-L1-R-SSR','RF-RR', 'dCRT (L1)', 'dCRT (RF)')),
    pvalue=CI_test_result$pvalues[ind]
  )
  temp_table=temp_table%>%
    mutate(setT=char_setT,pvalue=sprintf("%.3f", pvalue))
  return(temp_table)
}


load('seed2023_DIA_2023-11-12.RData')
CI_table_weak=CI_from_data_to_frame(CI_test_result_T1, char_setT="$\\mc{T}_1$")
CI_table_strong=CI_from_data_to_frame(CI_test_result_T2, char_setT="$\\mc{T}_2$")

CI_table=rbind(CI_table_weak, CI_table_strong)

CI_table$setT=as.factor(CI_table$setT)


tab=tabular( setT ~ Method*pvalue*c ,CI_table)


colLabels(tab)=colLabels(tab)[1:2,]
attr(rowLabels(tab),'dimnames')[[2]]='\\mc{T}'
latex.tabular(tab,file='Stock_pvalue.tex',sepC = 4, sepR = 2)


###########################  noise  power curve      ########################


library(reshape)
library(ggplot2)




load('Stock-noise-DIA-significant.RData')

noisy_pvalues_table=simplify2array(lapply(noisy_results, function(x)
  sapply(x, `[[`, 'pvalues')))

sig.lv=0.05

noise_power_table=apply(noisy_pvalues_table,c(3,1), function(v){mean(v<sig.lv) } ) 
noise_power_table=data.frame(noise_power_table)

noise_se_table=apply(noisy_pvalues_table,c(3,1), function(v){Fn_SE(v<sig.lv) } )
noise_se_table=data.frame(noise_se_table)

noise_se_table$noise_leve = noise_power_table$noise_leve = 
  as.numeric(rownames(noise_power_table))

colnames(noise_se_table)=colnames(noise_power_table)=
  c('LM-L1-R-SSR', 'RF-RR', 'dCRT (L1)', 'dCRT (RF)','noise_level')


df_melt=melt(noise_power_table,  id.vars="noise_level",  variable_name = "Test")
df_melt$Power=df_melt$value;df_melt$value=NULL

SE_melt <- melt(noise_se_table, id.vars = "noise_level", variable_name = "Test")
SE_melt$SE=SE_melt$value;SE_melt$value=NULL



# Merge the data and SE data frames
merged_df <- merge(df_melt, SE_melt, by = c("noise_level", "Test"))

merged_df$Test=factor(merged_df$Test, levels=c('LM-L1-R-SSR', 'RF-RR', 'dCRT (L1)', 'dCRT (RF)'))
merged_df=arrange(merged_df, Test)


num_methods=4
distinct_colors <- scales::hue_pal()(num_methods)
distinct_linetypes <- rep(c("solid", "dashed", "dotted", "twodash", "longdash", "dotdash"), length.out = num_methods)
distinct_shapes <- rep( c(1,2,3,4,5,6,8,15,16,17,18), length.out= num_methods)
distinct_linewide <- rep(c(0.6, 0.6, 1.2, 0.6, 0.6, 0.6), length.out = num_methods)

# Create ggplot visualization
gp = ggplot(merged_df, aes(x = noise_level, y = Power, color = Test, shape = Test, 
                           linetype = Test, linewidth = Test)) +
  geom_line() +
  scale_linetype_manual(values = distinct_linetypes) +
  scale_linewidth_manual(values = distinct_linewide) +
  scale_color_manual(values = distinct_colors) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1.2, 'cm'), 
        legend.key.width = unit(1.5, 'cm'),
        text = element_text(size = 18), # Adjusts overall text size; you can change the value as needed
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.title.x = element_text(size = 18), # Adjusts x-axis label size
        axis.title.y = element_text(size = 18), # Adjusts y-axis label size
        plot.title = element_text(size = 18) # Adjusts plot title size
  ) +
  geom_point(size = 3) +
  labs(x = expression(sigma),
       y = "Rejection rate")  +
  xlim(0.025, 0.25) + 
  ylim(0, 1)   + 
  scale_shape_manual(values = distinct_shapes) +
  scale_fill_manual(values = distinct_colors)  

print(gp)


new_gp <-  gp +
  geom_hline(aes(yintercept = 0.05, size = "Level = .05"), color = "black", linetype = "dotted") +
  scale_size_manual(name = "Parameters", values = c(0.5, "Level = .05" = 0.5)) +
  guides(size = guide_legend(title = "\n "))


# Print the plot
print(new_gp)

ggsave(filename = paste0("Stock-Noise.pdf"), plot = new_gp, device = "pdf", width = 8, height = 6)

