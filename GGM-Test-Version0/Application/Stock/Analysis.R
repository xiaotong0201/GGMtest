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



# Define X  ------

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

## Data transformation ---- 
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

## save(industry,d30_name,d30_ind,start_date, end_date, X,X_array,average_every_five, file='Pre-X.RData')

# Define Y-----



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


# Analysis begin ------

n = nrow(X)
p = ncol(X)

## Define and draw GICS----

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

plot(
  g,
  vertex.shape = "circle" ,
  vertex.size = 5,
  vertex.label.cex = 0.5,
  edge.color = "black",
  label.col='white',
  edge.width = 0.5
)

# Adding a legend
legend("bottomleft", legend = levels(industry), 
       col = cols, pch = 21, pt.cex = 2, cex = 1, bty = "n")



sort(colSums(G0), decreasing = T)
# max=18

### Estimation by Glasso -----
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
# 0.15 max=30 accept
# median 19 

sum(G1) / p / (p - 1) #18.4%


g <- graph.adjacency(G1 , mode = "undirected", diag = FALSE)

V(g)$color=cols[industry]
V(g)$label=stocks_name

pdf('stocks_G.pdf', width=8, height=6)

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

dev.off()


##   GoF test -----

source('../US/FunctionsUS.R')
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

# CI test----- 



wd_temp = getwd()
setwd('../../CRT')
source('CRT-Stat.R')
source('CRT-Computing.R')
source('Run-Testing.R')
setwd(wd_temp)

setting = 'a_st'



List_Stat_Comb=c(
  "Stat_Comb_Lasso_Res_Linear_Dev", 
  "Stat_Comb_RF_R_Res_RF"    
  , "Stat_Comb_RF_R_Distill_RF"
)

# Define experiment's X and Y ####
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
# save(G1, setT1, setT2, mat_X, Y, file="DataForExp.RData")




