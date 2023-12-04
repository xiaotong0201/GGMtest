# US Precipitation Analysis Script
# ---------------------------------------------------------

# Clearing the workspace
rm(list = ls()) # Removes all objects from the workspace

source('FunctionsUS.R')



# Preliminary

# Reading precipitation data
rain <- readRDS(file = "rain_state_year.rds")

# transform
x <- rain
# Box-Cox Transformation for each column
for (i in 1:ncol(rain)) {
  adjusted_column <- rain[, i] 
  adjusted_column <- adjusted_column - min(adjusted_column) + (1e-2)*(max(adjusted_column) - min(adjusted_column))
  boxcox_model <- MASS::boxcox(lm((adjusted_column) ~ 1),plotit = F)
  optimal_lambda <- boxcox_model$x[which.max(boxcox_model$y)]
  transformed_data <- (adjusted_column ^ optimal_lambda - 1) / optimal_lambda
  x[, i] <- transformed_data
}

# Shapiro-Wilk Test for Normality on each column
shapiro_p_values <- apply(x,2, function(column) shapiro.test(column)$p.value)
head(sort(shapiro_p_values))

# Assuming `data` is your dataset
ks_p_values <- apply(x,2, function(data){ks.test(data, "pnorm", mean(data), sd(data))$p.value})
head(sort(ks_p_values))




# Reading adjacency graph data
graph0 <- readRDS('us_stat_adj.rds')

# Removing specific columns and rows
ind <- c(which(colnames(graph0) == 'District of Columbia'), which(colSums(graph0) == 0))
G0 <- graph0[-ind, -ind]

# Aligning data with the graph structure
X <- x[, match(colnames(G0), colnames(x))]


p=ncol(G0)
max.union.degree=0
for(i in 1:p){for(j in 1:p){
  if(G0[i,j]==0&i!=j){
    max.union.degree=max(max.union.degree, length(union(which(G0[i,]==1), which(G0[j,]==1))))
  }
}}
print(max.union.degree)

# Defining list of CSS statistics
List_CSS_stat <- c("PRC_SS", "ERC_Z_SS", "F_max", "F_sum", "glr_glasso")
CSS_param = list(M = 100, L = 3)



# Running all tests
set.seed(2023)
result_full <- run_all_gof_tests(X, G0, List_CSS_stat = List_CSS_stat, 
                                 CSS_param=CSS_param,
                                 CSS = NULL, dir_prefix = 'ExperimentResults/',
                                 bSave=T, bReturnCSS=T
                                 )

save(X, shapiro_p_values, G0,List_CSS_stat, CSS_param, result_full, 
     file=paste0('ExperimentResults/Full-test','.RData'))


summarize_vector=function(result_list){
vec_result=c(result_list$VV_pvalue, result_list$DP_pvalue, 
             unlist(result_list$CSS_test_pvalues))
vec_result=vec_result[-5]
names(vec_result)=c('VV','DP','PRC','ERC', 'F$_{\\Sigma}$',"GLR-${\\ell_1}$" )
return(vec_result)
}



## timing
## 16s for CSS 400
## 20 for PRC, 5.5 for ERC, 12 for F, 4.7 for glr
## 60s one run


library(parallel)

###########  Simulation Type-I check -------


library(glasso)
p <- ncol(X)
N <- nrow(X)

# Initialize rho
rho <- matrix(0, nrow = p, ncol = p)
rho[G0 == 0] <- 1e8 * N * p
diag(rho) <- 0

# Glasso
GLS <- glasso(var(X), rho, penalize.diagonal = FALSE)
estimated_covariance <- GLS$w
estimated_mean=colMeans(X)

run_simulate=function(epo){
  source('FunctionsUS.R')
  set.seed(epo)
  
  simu_X = mvrnorm(N, estimated_mean, estimated_covariance)
  return(run_all_gof_tests(simu_X, G0, List_CSS_stat,  CSS_param = list(M = 100, L = 3),
                           CSS=NULL, bSave=F, bReturnCSS=F,'ExperimentResults/'))
}

subsampling_results=list()
num_repeats=40

# Detect the number of cores
no_cores <- detectCores() - 1  # leave one core free for system processes

# Create a cluster
cl <- makeCluster(no_cores)


  cat("Processing simulation:\n")
  
  clusterExport(cl, c('run_simulate',"num_repeats",
                      "estimated_covariance", 'estimated_mean',
                      'G0','N','p','List_CSS_stat')
  )
  
simulation_results <- parLapply(cl, 1:num_repeats , run_simulate)
  

save(simulation_results, num_repeats,estimated_covariance, file='ExperimentResults/US-null-simulation.RData')


load('ExperimentResults/US-null-simulation.RData')
error_rejection=sig.lv > t(sapply(simulation_results, summarize_vector))
error=colMeans(error_rejection)
se_error=apply(error_rejection, 2, function(v)sd(v)/sqrt(length(v)))


simu_long_table=data.frame(Method=as.factor(names(error)),
                      Power=paste(sprintf("%.3f", error), sprintf("(%.3f)", se_error)))

pvalues=summarize_vector(result_full)
pv_long_table=data.frame(Method=as.factor(names(pvalues)),pvalues=sprintf("%.3f",  pvalues))

merged_table=merge(simu_long_table,pv_long_table,by='Method')
merged_table=merged_table[match(names(pvalues), merged_table$Method  ), ]
merged_table$Method = factor( merged_table$Method , levels= names(pvalues))

tab=tabular( (pvalues+Power)  ~ Method*c ,merged_table)
colLabels(tab)=colLabels(tab)[1:2,]
rowLabels(tab)[1,1]='P-Value'
rowLabels(tab)[2,1]='Type-I error'

latex.tabular(tab,file='US.tex',sepC = 6, sepR = 1)


tab=tabular( Power  ~ Method*c ,merged_table)
colLabels(tab)=colLabels(tab)[1:2,]
rowLabels(tab)[1,1]='Type-I error'
latex.tabular(tab,file='US_type_i_err.tex',sepC = 6, sepR = 1)



tab=tabular( pvalues  ~ Method*c ,merged_table)
colLabels(tab)=colLabels(tab)[1:2,]
rowLabels(tab)[1,1]='P-Value'
latex.tabular(tab,file='US_pvalues.tex',sepC = 6, sepR = 1)



########## subsampling -------
## 50 minutes for one subsample size 

# Define sample sizes and number of repeats
sample_sizes <- c(18, 21, 24, 27, 30) # seq(18, 30, by = 2)
num_repeats <- 400


run_subsample=function(epo){
  source('FunctionsUS.R')
  set.seed(epo)
  N=nrow(X)
  subsample <- X[sample(N, sn), ]
  return(run_all_gof_tests(subsample, G0, List_CSS_stat,  CSS_param = list(M = 100, L = 3),
                           CSS=NULL, bSave=F, bReturnCSS=F,'ExperimentResults/'))
}

subsampling_results=list()

# Detect the number of cores
no_cores <- detectCores() - 1  # leave one core free for system processes

# Create a cluster
cl <- makeCluster(no_cores)

for (sn in sample_sizes) {
  cat("Processing sample size:", sn, "\n")
  
  clusterExport(cl, c("run_subsample", "num_repeats",
                      'G0','X','sn','List_CSS_stat')
  )
  
  new_mc_results <- parLapply(cl, 1:num_repeats , run_subsample)
  subsampling_results[[as.character(sn)]] <- new_mc_results
}
stopCluster(cl)

save(subsampling_results, num_repeats,sample_sizes, file='ExperimentResults/US-subsampling.RData')


