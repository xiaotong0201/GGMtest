library(dplyr)
library(tables)
source("~/Library/CloudStorage/GoogleDrive-statdongminghuang@gmail.com/My Drive/Codes in R/ArrayToTable-Package.R")


########  Type 1 error ######


pvalues=t(sapply(null_result, `[[`, "pvalues"))
pvalues_lrt=pvalues[,1]
# Create the histogram
bw_histogram <- ggplot(data.frame(x = pvalues_lrt), aes(x)) +
  geom_histogram(bins = 20, fill = "grey", color = "white") + # Dark bars with light borders
  theme_minimal() + # Minimalist theme
  theme(text = element_text(size = 20), # Adjust text size for readability
        axis.title = element_text(size = 18), # Slightly larger axis titles
        plot.title = element_text(size = 20, hjust = 0.5)) + # Centered title
  labs(title = "Histogram",
       x = "Simulated p-value",
       y = "Frequency")
print(bw_histogram)


ggsave("P-LRT_histogram.pdf", bw_histogram, width = 8, height = 6)


Fn_SE=function(x){xx=x[!is.na(x)]; sd(xx)/sqrt(length(xx))}

null_from_data_to_frame=function(null_result){
  
  
  
  null_pvalue_table = t(sapply(null_result,  '[[', "pvalues"))
  power=apply(null_pvalue_table, 2, function(x)mean(x<0.05, na.rm=T))
  se=apply(null_pvalue_table,2,Fn_SE)
  temp_table=data.frame(
    Method=factor(c('$G$-CRT (RF)','LRT', 'dCRT (L1)', 'dCRT (RF)'),
                  levels = c('LRT', '$G$-CRT (RF)','dCRT (L1)', 'dCRT (RF)') ),
    power=power, 
    se=se
  )
  temp_table=temp_table%>%
    mutate(Error=sprintf("%.3f (%.3f)", power, se))
  return(temp_table)
}

load('ERplus-Type-I.RData')

type1_table=null_from_data_to_frame(null_result)


tab=tabular(  ~ Method*Error*c ,type1_table)


colLabels(tab)=colLabels(tab)[1:2,]
latex.tabular(tab,file='BCR_Type-I-error.tex',sepC = 4, sepR = 1)



