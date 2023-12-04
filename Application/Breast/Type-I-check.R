library(parallel)
### type i error check
source('~/Library/CloudStorage/GoogleDrive-statdongminghuang@gmail.com/My Drive/Codes in R/Template_foreach_evendistributed.R')


# load("DataForExp.RData")
load('ERplus-CI-testing-seed2023.RData')

prefunc = function() {
  setwd('../../CRT')
  source('../GGMTesting.R')
  source('../CSS-GGM.R')
  source('CRT-Stat.R')
  source('CRT-Computing.R')
  source('Run-Testing.R')
  library(VGAM)
  
}


############# baseline



null_result = list()
null_func = function(i) {
  set.seed(i)
  
  
  sim_X = rmvnormal(nrow(mat_X), mu = colMeans(mat_X), sigma = cov(mat_X))

  fit_value_null=predict(lasso_null,newx = sim_X[, setV], type='response', s=cv_lasso_null$lambda.min)
  Y_sim=t(apply(fit_value_null, 1, function(v)rmultinom(1,1,v)))
  Y_sim_factor=apply(Y_sim==1,1,which)
  
  CSS_param = list(M = 400, L = 3)
  CSS_sample = exchangeable_sampling(
    sim_X,
    G1,
    I = setT,
    bReduced = T,
    M = CSS_param$M,
    L = CSS_param$L
  )

CI_result = Run_all_test(
  sim_X,
  Y_sim_factor, 
  setT,
  CSS_sample,
  setting = 'a_bc',
  List_Stat_Comb = "Stat_Comb_RF_C_IS", 
  IncludeBaseline = T, 
  )
  
  return(CI_result)
  
}


set.seed(2023)
setV=setdiff(1:ncol(mat_X), setT)


cv_lasso_null=cv.glmnet(mat_X[, setV], sub_Y,family='multinomial',nfolds = 5)
lasso_null=glmnet(mat_X[, setV], sub_Y,family='multinomial')


new_result = run_parallel(
  replication_start = 1,
  replication_end = 400,
  num_cores = 7,
  list_export = c(
    'cv_lasso_null','lasso_null',
    'mat_X',
    'G1',    'setT', 'setV',
    'run_and_time',
    'baseline_method_calls',
    'List_Stat_Comb'
  ),
  func_name = 'null_func',
  prefunc_name = 'prefunc'
)

null_result=c(null_result, new_result)

# save(null_result, file='ERplus-Type-I.RData')

load('ERplus-Type-I.RData')


pvalues=t(sapply(null_result, `[[`, "pvalues"))

apply(pvalues, 2, function(x)mean(x<0.05, na.rm=T))


apply(pvalues, 2, function(x)mean(x<0.05, na.rm=T))
apply(pvalues<0.05, 2,Fn_SE)

pvalues_lrt=pvalues[,2]

hist(pvalues_lrt)


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

# Print the plot
print(bw_histogram)

# Optionally, save the plot
ggsave("P-LRT_histogram.pdf", bw_histogram, width = 8, height = 6)
