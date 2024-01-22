# Loading required libraries
source('Function.R')
source('FunctionsUS.R')
library(tables)
library(reshape)
library(plyr)
library(dplyr)
library(stringr)
library(Matrix)
library(MASS)
library(mvtnorm)
library(glasso)
library(parallel)
library(Rfast)
library(glmnet)
library(rags2ridges)
library(tidyquant)

##########################################################################
        ###### Average Daily Precipitation in the United States ######
##########################################################################

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

################## Create Graph ################## 

# Mapping between two-letter abbreviations and full state names
state_mapping <- c(
  "AL" = "Alabama", "AK" = "Alaska", "AZ" = "Arizona", "AR" = "Arkansas", "CA" = "California",
  "CO" = "Colorado", "CT" = "Connecticut", "DE" = "Delaware", "DC" = "District of Columbia",
  "FL" = "Florida", "GA" = "Georgia", "HI" = "Hawaii", "ID" = "Idaho", "IL" = "Illinois",
  "IN" = "Indiana", "IA" = "Iowa", "KS" = "Kansas", "KY" = "Kentucky", "LA" = "Louisiana",
  "ME" = "Maine", "MD" = "Maryland", "MA" = "Massachusetts", "MI" = "Michigan", "MN" = "Minnesota",
  "MS" = "Mississippi", "MO" = "Missouri", "MT" = "Montana", "NE" = "Nebraska", "NV" = "Nevada",
  "NH" = "New Hampshire", "NJ" = "New Jersey", "NM" = "New Mexico", "NY" = "New York",
  "NC" = "North Carolina", "ND" = "North Dakota", "OH" = "Ohio", "OK" = "Oklahoma", "OR" = "Oregon",
  "PA" = "Pennsylvania", "RI" = "Rhode Island", "SC" = "South Carolina", "SD" = "South Dakota",
  "TN" = "Tennessee", "TX" = "Texas", "UT" = "Utah", "VT" = "Vermont", "VA" = "Virginia",
  "WA" = "Washington", "WV" = "West Virginia", "WI" = "Wisconsin", "WY" = "Wyoming"
)

# Read the file
file_path <- 'state_neighbors.txt'
lines <- readLines(file_path)

# Initialize an empty adjacency matrix
n_states <- length(state_mapping)
adj_matrix <- matrix(0, nrow = n_states, ncol = n_states, 
                     dimnames = list((state_mapping), (state_mapping)))

# Populate the adjacency matrix
for(line in lines) {
  tokens <- unlist(strsplit(line, " "))
  state <- state_mapping[tokens[1]]
  neighbors <- state_mapping[tokens[-1]]
  
  for(neighbor in neighbors) {
    adj_matrix[state, neighbor] <- 1
    adj_matrix[neighbor, state] <- 1  # Assuming undirected graph
  }
}

# Display the adjacency matrix
#print(adj_matrix)
graph0=adj_matrix
saveRDS(graph0,file = 'us_stat_adj.rds')

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

################## Simulation Type-I check ################## 

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
  library(MASS)
  simu_X = mvrnorm(N, estimated_mean, estimated_covariance)
  
  
  return(run_all_gof_tests(simu_X, G0, List_CSS_stat,  CSS_param = list(M = 100, L = 3),
                           CSS=NULL, dir_prefix ='ExperimentResults/',bSave=T, bReturnCSS=T))
}

subsampling_results=list()
num_repeats=400

# Detect the number of cores
no_cores <- detectCores() - 1  # leave one core free for system processes

# Create a cluster
cl <- makeCluster(no_cores)


#cat("Processing simulation:\n")

clusterExport(cl, c('run_simulate',"num_repeats",
                    "estimated_covariance", 'estimated_mean',
                    'G0','N','p','List_CSS_stat')
)

simu_X=mvrnorm(N, estimated_mean, estimated_covariance)
run_all_gof_tests(simu_X, G0, List_CSS_stat,  CSS_param = list(M = 100, L = 3),
                  CSS=NULL, dir_prefix ='ExperimentResults/',bSave=T, bReturnCSS=T)

simulation_results <- parLapply(cl, 1:num_repeats , run_simulate)

save(simulation_results, num_repeats,estimated_covariance, file='ExperimentResults/US-null-simulation.RData')

load('ExperimentResults/US-null-simulation.RData')
sig.lv=0.05
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
print(tab)

################## Subsampling ##################

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

##########################################################################
               ###### Dependence of Fund Return ######
##########################################################################

# Load necessary libraries
library(tidyquant)

List_stocks=read.csv('list_stocks.csv')
stocks_name=List_stocks$Symbol
industry=as.factor(List_stocks$Sector)


# Define the list of Dow 30 stock names
d30_name <- c(
  "AAPL", "AMGN", "AXP", "BA", "CAT","CRM", "CSCO","CVX", "DIS","DOW",   "GS",
  "HD", "HON", "IBM", "INTC", "JNJ", "JPM", "KO","MCD", "MMM","MRK", "MSFT", "NKE", "PG","TRV",
  "UNH",  "V","VZ", "WBA", "WMT"
)

# Match Dow 30 stock names with S&P 100 stock names
d30_ind <- match(d30_name, stocks_name)
names(d30_ind)=d30_name

summary(industry)
summary(industry[d30_ind])

################## Define X and Y ##################

# Define X

if (!file.exists('stocks-2020-09-01-20231030.RData')) {
  # Download data for S&P 100 stocks
  X_list <- list()
  for (i in 1:length(stocks_name)) {
    stock_name <- stocks_name[i]
    stock_data <-
      tq_get(stock_name, from = '2020-09-01', to = '2023-10-30')
    X_list[[stock_name]] <- stock_data
  }
  save(X_list, stocks_name, file = 'stocks-2020-09-01-20231030.RData')
} else{
  load(file = 'stocks-2020-09-01-20231030.RData')
}

# Define the start and end dates for the selection
start_date <- as.Date('2020-09-01')
end_date <- as.Date('2022-06-30')

# Extract and process the stock data
X_array <- sapply(X_list, function(v) {
  selected <- v[which((v$date >= start_date) & (v$date <= end_date)), ]
  u <- c(selected$adjusted)
  names(u) <- selected$date
  return(u)
})


# Calculate log returns and simple returns
log_return <- apply(X_array, 2, function(v) {
  diff(log(v))
})
return <- exp(log_return) - 1

# Plot histograms and a scatter plot
hist(return)
hist(log_return)
plot(log_return ~ return)  # on a line

## Data transformation 
# Function to calculate the average of every 5 entries for a column
average_every_five <- function(column) {
  group <- (seq_along(column) - 1) %/% 5 + 1
  tapply(column, group, mean)
}

# Apply the function to each column of the data frame
ave_return <- apply(return, 2, average_every_five)

# Transform the average returns and tune parameters for each column
source('sym_piece_power_log.R')
transformed_AveReturn <- ave_return
tuned_Parameters <- list()

for (i in 1:ncol(log_return)) {
  column_Adjusted <- ave_return[, i]
  column_Adjusted <- column_Adjusted - median(column_Adjusted)
  best_Parameters <- find_best_parameters(
    column_Adjusted,
    range_lambda = c(0.1, 2),
    range_gamma = quantile(abs(column_Adjusted), c(0.8, 1)) 
  )
  tuned_Parameters[[i]] <- best_Parameters
  transformed_Column <- transform_piece_data(column_Adjusted, best_Parameters[1], cutoff = best_Parameters[2])
  transformed_AveReturn[, i] <- transformed_Column
}

# Apply Shapiro-Wilk Test for normality on each column of the transformed data
sym_shapiro_p_values <-
  apply(transformed_AveReturn, 2, function(column)
    shapiro.test(column)$p.value)
head(sort(sym_shapiro_p_values))

# Box-Cox Transformation for each column in average returns
boxCox_TransformedData <- ave_return
for (i in 1:ncol(ave_return)) {
  column_Adjusted <- ave_return[, i]
  column_Adjusted <-
    column_Adjusted - min(column_Adjusted) + 1e-2 * (max(column_Adjusted) - min(column_Adjusted))
  boxcox_Model <-
    MASS::boxcox(lm(column_Adjusted ~ 1), plotit = FALSE)
  optimal_Lambda <- boxcox_Model$x[which.max(boxcox_Model$y)]
  transformed_Column <-
    (column_Adjusted ^ optimal_Lambda - 1) / optimal_Lambda
  boxCox_TransformedData[, i] <- transformed_Column
}
bc_shapiro_p_values <-
  apply(boxCox_TransformedData, 2, function(column)
    shapiro.test(column)$p.value)
head(sort(bc_shapiro_p_values))


plot(log(sort(sym_shapiro_p_values)))
points(log(sort(bc_shapiro_p_values)), add = T)
abline(h = log(0.1), col = 'red')

plot(
  sym_shapiro_p_values,
  bc_shapiro_p_values,
  xlim = c(0, 0.2),
  ylim = c(0, 0.2)
)

# Identify columns failing the Shapiro-Wilk test for symmetry
ind_sym_fail = which(sym_shapiro_p_values < 0.05)
if( min( bc_shapiro_p_values[ind_sym_fail] ) < 0.05){
  stop('Both transformation failed')
}

# Combine transformed data, substituting columns that failed the symmetry test
final_TransformedData <- transformed_AveReturn
final_TransformedData[, ind_sym_fail] <- boxCox_TransformedData[, ind_sym_fail]
X <- scale(final_TransformedData)
hist(X)
qqnorm(X)
qqline(X, col = 'red')

# Define Y

DIA = getSymbols(
  'DIA',
  from = start_date,
  to = end_date,
  warnings = FALSE,
  src = "yahoo",
  auto.assign = FALSE
)

data_y = DIA[, 6]

dim(data_y)
return_ndji = exp(diff(log(data_y)) [-1]) - 1
Y = c(average_every_five(return_ndji)) *100


length(Y) == nrow(X)


# Analysis begin

n = nrow(X)
p = ncol(X)

## Define and draw GICS

G0 = matrix(0, p, p)


for (i in 1:p) {
  for (j in 1:p) {
    if (industry[i] == industry[j]) {
      G0[i, j] = 1
    }
  }
}
library(igraph)

g <- graph.adjacency(G0, mode = "undirected", diag = FALSE)
cols <-
  c(
    "tomato1",
    "orange",
    "gold1",
    "darkolivegreen1",
    "lightblue",
    "slateblue2",
    "mediumorchid2",
    "green",
    "lightpink1",
    "royalblue4",
    "cyan"
  )

sectors = unique(industry)

V(g)$color=cols[industry]
V(g)$label=stocks_name

sort(colSums(G0), decreasing = T)


################## Estimation by Glasso ##################

library(glasso)


rho0 = 1 - G0
diag(rho0) = 0
rho = rho0 *  0.15 


GLS <- glasso(var(X), rho = rho)
inv.sigma0_gls <- GLS$wi #inverse covariance matrix

G1 = (inv.sigma0_gls != 0)
diag(G1) = 0

median(colSums(G1))
table(colSums(G1))
# 0.15 max=29 accept
# median 19 

sum(G1) / p / (p - 1) #18.4%


g <- graph.adjacency(G1 , mode = "undirected", diag = FALSE)

V(g)$color=cols[industry]
V(g)$label=stocks_name



# myigraphplot(
set.seed(1)
plot(
  g,
  vertex.shape = "circle" ,
  vertex.size = 10,
  vertex.label.color='black',
  vertex.label.cex = 0.35,
  edge.color = "black",
  edge.width = 0.8
)

# Adding a legend
legend("bottomright", legend = levels(industry), 
       col = cols, 
       pch = 19, 
       pt.cex = 1,   # Adjust point size
       cex = 0.5,    # Adjust text size
       bty = "o",    # Type of box. Use "n" for no box
       ncol=2,
       # nrow=4,
       lwd = 0)  

################## GoF Test ##################

set.seed(2023)
List_CSS_stat <- c("F_sum")
CSS_param = list(M = 1000, L = 3)

GoF_result <-
  run_all_gof_tests(
    X,
    G1,
    List_CSS_stat = List_CSS_stat,
    CSS_param = CSS_param,
    CSS = NULL,
    dir_prefix = 'ExperimentResults/',
    bSave = F,
    bReturnCSS = T
  )
GoF_result[1:3]

save(X, G1,CSS_param,GoF_result, file='G1GoF.Rdata'  )

################## Type I Error Check ##################

setting = 'a_st'

List_Stat_Comb=c(
  "Stat_Comb_Lasso_Res_Linear_Dev", 
  "Stat_Comb_RF_R_Res_RF"    
  , "Stat_Comb_RF_R_Distill_RF"
)

# Define experiment's X and Y 
mat_X = X 
dimnames(mat_X) = NULL

setT1 = match( c( "CAT", "DOW",'HON','MMM','PG') , stocks_name) 
setT2 = c(8,  28,  70, 94)  

############# baseline

baseline_method_calls = list(
  DebiasLasso = quote(DL_hdi(x, y, setT, model_type = "gaussian")),
  dCRT_Lasso = quote(dCRT(x, y, setT, model_type = "Gaussian_lasso")),
  dCRT_RF = quote(dCRT(x, y, setT, model_type = "RF"))
)

List_Stat_Comb = c(
  "Stat_Comb_Lasso_Res_Linear_Dev",
  "Stat_Comb_RF_R_Res_RF"
)

setT=setT2

setV=setdiff(1:ncol(mat_X), setT)
cv_lasso_null=cv.glmnet(mat_X[, setV], Y)
fit_value_null=predict(cv_lasso_null,newx = mat_X[, setV], type='response', s=cv_lasso_null$lambda.min)
err_sd_null=sd(Y-fit_value_null)
coef_null = predict(cv_lasso_null,type='coef', s=cv_lasso_null$lambda.min)

null_result = list()
for (i in 1:400){
  set.seed(i)
  
  # Generate simulated data
  sim_sample = exchangeable_sampling(
    mat_X,
    G1, 
    bReduced = F,
    M = 1,
    L = 20
  )
  sim_X=sim_sample$X_copies[[1]] 
  Y_sim =coef_null[1]+sim_X[,-setT]%*% c(coef_null[-1]) + err_sd_null * rnorm(nrow(mat_X))
  
  CSS_param = list(M = 400, L = 3)
  CSS_sample = exchangeable_sampling(
    sim_X,
    G1,
    I = setT,
    bReduced = T,
    M = CSS_param$M,
    L = CSS_param$L
  )
  
  CI_test_result = Run_all_test(
    sim_X,
    Y_sim,
    setT,
    CSS_sample,
    setting = 'a_st',
    List_Stat_Comb = List_Stat_Comb,
    IncludeBaseline = T
  )
  
  null_result[[i]]=CI_test_result
  
}

pvalues=t(sapply(null_result, `[[`, "pvalues"))
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

type1_table_strong=null_from_data_to_frame(null_result, char_setT="$\\mc{T}_2$")

################## CIT Test ##################

setting = 'a_st'

List_Stat_Comb=c(
  "Stat_Comb_Lasso_Res_Linear_Dev", 
  "Stat_Comb_RF_R_Res_RF"    
  , "Stat_Comb_RF_R_Distill_RF"
)

# Define experiment's X and Y 
mat_X = X 
dimnames(mat_X) = NULL

setT1 = match( c( "CAT", "DOW",'HON','MMM','PG') , stocks_name) 
setT2 = c(8,  28,  70, 94)  


exp_seed1=202311
set.seed(exp_seed1)

CSS_param = list(M = 1000, L = 3)
CSS_sample_T1 = exchangeable_sampling(
  mat_X,
  G1,
  I = setT1,
  bReduced = T,
  M = CSS_param$M,
  L = CSS_param$L
)

CI_test_result_T1 = Run_all_test(
  mat_X,
  Y,  
  setT1,
  CSS_sample_T1,
  setting = 'a_st',
  List_Stat_Comb =  List_Stat_Comb,  
  IncludeBaseline = T
)


exp_seed2=20231112
set.seed(exp_seed2)

CSS_param = list(M = 1000, L = 3)
CSS_sample_T2 = exchangeable_sampling(
  mat_X,
  G1,
  I = setT2,
  bReduced = T,
  M = CSS_param$M,
  L = CSS_param$L
)

CI_test_result_T2 = Run_all_test(
  mat_X,
  Y,  
  setT2,
  CSS_sample_T2,
  setting = 'a_st',
  List_Stat_Comb =  List_Stat_Comb,  
  IncludeBaseline = T
)

if('Saving'==T){
  save(G1, mat_X,Y, setT1,exp_seed1, CSS_sample_T1,CI_test_result_T1, setT2,exp_seed2, CSS_sample_T2, CI_test_result_T2,  
       file='seed2023_DIA_2023-11-12.RData')
}

##########################################################################
                   ###### Breast Cancer Relapse ######
##########################################################################

# Y part 
# Load the CSV data
table_path = 'GSE2034-PI.csv'
PI = read.csv(table_path)

# Display the first few rows of the loaded data
head(PI)

# subclass
ER=factor(PI$ER.Status,levels=c("ER+","ER-"))

# Survival time
ST=PI$time.to.relapse.or.last.follow.up..months.
NotCensor=PI$relapse..1.True.

# Create a binary variable Y.brain based on conditions
Y.brain = rep(0, nrow(PI))
Y.brain[which(PI[, 4] <= 48)] = 1
Y.brain[which(PI[, 4] <= 48 & PI[, 7])] = 2
Y.brain = as.factor(Y.brain)

# Display a summary of Y.brain
summary(Y.brain)
Y=Y.brain

# X part 
# Load gene data
Data = read.csv('GSE2034_series_matrix.txt', skip = 26, sep = '\t')

GeneData = t(Data[-(1:28),])
GeneName = GeneData[1,]
GeneData = apply(GeneData[-1,],MARGIN = c(1, 2), FUN = as.numeric)

ImpVarName = read.csv('ImpVar.csv',header = F)$V1
ImpVarInd =  which(GeneName %in% ImpVarName)  #match(ImpVarName,GeneName) #
X_gene = GeneData[
  match(PI$GEO.asscession.number,Data[1,-1])
  ,ImpVarInd]


# Display the dimensions of X
dim(X_gene)

# Create a histogram for the first column of X
hist(X_gene[,2])


library(bestNormalize)
X_norm=apply(X_gene, 2, function(v)orderNorm(v)$x.t)
order_shapiro_p_values <-
  apply(X_norm, 2, function(column)
    shapiro.test(column)$p.value)
sort(order_shapiro_p_values, index=T)


### check if the mariginal is significant
library(survival)
X=X_norm

for(j in 1:ncol(X)){
  linearmodel=glm((Y!=0)~X[,j],family=binomial)
  print(summary(linearmodel)$coef[2,4])
}

for(j in 1:ncol(X)){
  cox_model=coxph(Surv(time = ST,event = NotCensor) ~ X[,j])
  print(summary(cox_model)$coef[1,5])
}


library(glmnet)
# Goal: Define T  
all.rel.set=1:76

# May use ER
sum(ER=='ER+')

sub_X=X[ER=='ER+', 1:60 ]
sub_Y=Y[ER=='ER+']


# Method 1, use all stocks to run a Lasso selection
set.seed(2023)
lasso.cv=cv.glmnet(sub_X,sub_Y,family="multinomial",nfolds = 5)
lasso=glmnet(sub_X ,sub_Y,family="multinomial")
betahat_list=coef(lasso,s=lasso.cv$lambda.min)
betahat_mat=do.call(cbind,betahat_list)
betahat_nonzero=which((rowSums(betahat_mat^2)[-1])!=0)
strong.set=intersect(all.rel.set, betahat_nonzero)
weak.set=setdiff(all.rel.set,betahat_nonzero)


## Meaning: whether should use the best CV lambda or should 
##  include more predictors
setT_multin_weak=weak.set


plot(lasso)
betahat_ose_list=coef(lasso,s=lasso.cv$lambda.1se)
betahat_ose_mat=do.call(cbind,betahat_ose_list)
betahat_ose_nonzero=which((rowSums(betahat_ose_mat^2)[-1])!=0)

## Meaning: whether should use the best CV lambda or the 1se lambda
setT_multin_1se=setdiff(strong.set,betahat_ose_nonzero) 


GeneName[ImpVarInd[betahat_ose_nonzero]]
GeneName[ImpVarInd[setT_multin_1se]]

################## Type I Error ##################

mat_X=sub_X[, strong.set]
p=ncol(mat_X)
setT= match(setT_multin_1se , strong.set)
G1=matrix(1,p,p)
diag(G1)=0

set.seed(2023)
setV=setdiff(1:ncol(mat_X), setT)
cv_lasso_null=cv.glmnet(mat_X[, setV], sub_Y,family='multinomial',nfolds = 5)
lasso_null=glmnet(mat_X[, setV], sub_Y,family='multinomial')

null_result = list()
for (i in 1:400){
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
  
  null_result[[i]]=CI_result
  
}


pvalues=t(sapply(null_result, `[[`, "pvalues"))

###  likelihood ratio test
apply(pvalues, 2, function(x)mean(x<0.05, na.rm=T))
pvalues_lrt=pvalues[,2]
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

### Type-I error table
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

type1_table=null_from_data_to_frame(null_result)
print(type1_table)

################## CIT Test ##################

mat_X=sub_X[, strong.set]
p=ncol(mat_X)
setT= match(setT_multin_1se , strong.set)


#Chisq_Multi(mat_X, sub_Y, setT)


G1=matrix(1,p,p)
diag(G1)=0


exp_seed=2023
set.seed(exp_seed)

CSS_param = list(M = 1000, L = 3)

CSS_sample = exchangeable_sampling(
  mat_X,
  G1,
  I = setT,
  bReduced = T,
  M = CSS_param$M,
  L = CSS_param$L
)

## random forest ok
RF_result = Run_all_test(
  mat_X,
  sub_Y, 
  setT,
  CSS_sample,
  setting = 'a_bc',
  List_Stat_Comb = "Stat_Comb_RF_C_IS" ,
  IncludeBaseline = T, 
)
print(RF_result)

save(mat_X, sub_Y, setT, G1,CSS_sample, RF_result, file='ERplus-CI-testing-seed2023.RData' )

################## Noise Experiments ##################

load("DataForExp.RData")

sub_Y=Y 
setT=setT2

library(parallel)

noise_level= c(0.05, 0.10, 0.15, 0.20, 0.25)


num_repeats <- 400


run_noisy = function(epo) {
  source('FunctionsStock.R')
  set.seed(epo)
  
  
  exp_Y <- sub_Y + rnorm(length(sub_Y), sd = noise_sd)
  
  
  return(para_run_one(mat_X, exp_Y, setT, G1))
}

# run_noisy(1)   # test before running

#### Parallel begin
noisy_results = list()

# Detect the number of cores
no_cores <-
  detectCores()   # leave one core free for system processes

# Create a cluster
cl <- makeCluster(no_cores)

for (noise_sd in noise_level  ) {
  cat("Processing noise sd:", noise_sd, "\n")
  
  clusterExport(cl,
                c("run_noisy", "num_repeats",
                  'G1', 'mat_X', 'setT', 'sub_Y', 'noise_sd'))
  
  new_mc_results <- parLapply(cl, 1:num_repeats , run_noisy)
  noisy_results[[as.character(noise_sd)]] <- new_mc_results
}
stopCluster(cl)



save(noisy_results, num_repeats,noise_level, file = 'Stock-noise-DIA-significant.RData')


sum(unlist(lapply(noisy_results, function(x)sapply(x, `[[`, 'use.time')))) / 4


noisy_pvalues_table=simplify2array(lapply(noisy_results, function(x)
  sapply(x, `[[`, 'pvalues')))

sig.lv=0.05
apply(noisy_pvalues_table,c(3,1), function(v){mean(v<sig.lv) } )
