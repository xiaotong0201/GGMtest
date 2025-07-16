rm(list = ls())

library(Matrix)
library(MASS)
library(mvtnorm)
library(glasso)
library(Rfast)
library(glmnet)
library(rags2ridges)
source('../CSS-GGM.R')
source('../GGMTesting.R')
source('GoF-Statistics.R')
source('GoF-Settings.R')
source('OneExp_varyL.R')

FlagLocalTest = F

param_DenseER_template <- list(model = "ER",
                               N = c(50),
                               P = c(120),
                               modelParam = list(
                                 Q = c(0.5),
                                 Q0 = 0.25,
                                 s = c(0.01,0.02)
                               )
)

run_experiments_epo <- function(param_template, epo, M, L_values, List_CSS_stat = NULL,
                                dir_prefix = '') {
  final_results <- list()
  lFR = 0
  for (N in param_template$N) {
    for (P in param_template$P) {
      # Create a list to hold the possible combinations of model parameters
      combos <- expand.grid(param_template$modelParam)
      
      # Loop over each combination of model parameters
      for (i in 1:nrow(combos)) {
        modelParam = as.list(combos[i, ])
        names(modelParam) = names(param_template$modelParam)
        
        param_part <- list(model = param_template$model, 
                           N = N, P = P, modelParam = modelParam)
        
        # Loop over each value of L
        
          # Call GGM function and store the result
          final_results[[lFR + 1]] = GGM_varyL(param_part, 
                                               L_values = L_values,
                                         CSS_param = list(M = M),
                                         epo = epo, List_CSS_stat = List_CSS_stat,
                                         dir_prefix = dir_prefix)
          lFR = lFR + 1
        
      }
    }
  }
  return(final_results)
}

args = commandArgs(trailingOnly = TRUE)
ind.experiment = as.numeric(args[1])
epo = as.numeric(args[2])

# small test
if (FlagLocalTest) {
  seq_L <- c(1, 2, 3, 5, 10, 20)
  M <- 20
  
  # Define a small test parameter template
  test_param_template <- list(model = "ER",
                              N = c(50),
                              P = c(120),
                              modelParam = list(
                                Q = c(0.5),
                                Q0 = 0.25,
                                s = c(0.01,0.02)
                              )
  )
  
  param_template <- test_param_template
  ind.experiment <- 5
  epo=1
} else {
  param_template <- param_DenseER_template
  seq_L <- c(1, 2, 3, 5, 10, 20)
  M <- 100

  }

PlatForm = (Sys.info())["sysname"]
if (FlagLocalTest & (PlatForm != "Windows")) { Rprof("my_profile.out") }

path_prefix = paste0('../../../Programming/LocalExperiments/GGM/GoF/Experiments-', ind.experiment, '/')
if (!dir.exists(path_prefix)) {
  dir.create(path_prefix)
}

save(PlatForm, file = paste0(path_prefix, epo, '.RData'))

###############################################################################

output = run_experiments_epo(param_template, epo, M, L_values=seq_L, dir_prefix = path_prefix)
print(output)

save(output, param_template,
     file = paste0(path_prefix, epo, '.RData'))

###############################################################################

if (FlagLocalTest & (PlatForm != "Windows")) {
  Rprof(NULL)
  summary_data <- summaryRprof("my_profile.out")
  print(summary_data)
}
