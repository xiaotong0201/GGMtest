rm(list=ls())

library(Matrix)
library(MASS)
library(mvtnorm)
library(glasso)
library(Rfast)
library(glmnet)
library(rags2ridges)
source('../../CSS-GGM.R')
source('../../GGMTesting.R')
source('../GoF-Statistics.R')
source('MoreGraphs.R')


## Function to simulate the data
simulate_data <- function(Sigma, n) {
  X <- mvrnorm(n, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
  return(X)
}


run_experiment <- function(rep_id, n, p, model, dir_prefix = "Results", List_CSS_stat = NULL, CSS_param = list(M = 10, L = 5)) {
  if (is.null(List_CSS_stat)) {
    List_CSS_stat <- c(
      "PRC_SS", "PRC_SA", "PRC_Z_SS", "PRC_Z_SA",
      "ERC_SS", "ERC_SA", "ERC_Z_SS", "ERC_Z_SA",
      "F_max", "F_sum",
      "glr_glasso"
    )
  }
  
  # Create subdirectory for this replication
  model_dir <- file.path(dir_prefix, model)
  rep_dir <- file.path(model_dir, paste0("rep_", rep_id))
  if (!dir.exists(rep_dir)) dir.create(rep_dir, recursive = TRUE)
  
  # Setup graph and Sigma
  setup <- setup_graph_and_sigma(model, p)
  G0 <- setup$G0
  Sigma <- setup$Sigma
  
  if (min(eigen(Sigma)$values) <= 0) {
    cat("Error: Invalid covariance matrix.\n")
    return(NULL)
  }
  
  # Set seed for reproducibility
  set.seed(rep_id)
  
  # Simulate data
  x <- simulate_data(Sigma, n)
  saveRDS(x, file = file.path(rep_dir, "data.rds"))
  
  # Perform inference
  VV_pvalue <- VV_GoF(x, G0)
  DP_pvalue <- Bonf_GoF(x, G0)
  
  # CSS sampling
  CSS <- exchangeable_sampling(x, G0, M = CSS_param$M, L = CSS_param$L)
  saveRDS(CSS, file = file.path(rep_dir, "CSS.rds"))
  
  # CSS testing
  CSS_test_result <- lapply(List_CSS_stat, function(f) {
    CSSGoF(x, G0, CSS$X_copies, Fun_Stat = get(f))
  })
  saveRDS(CSS_test_result, file = file.path(rep_dir, "CSS_results.rds"))
  
  # Save RData for this replication
  save(VV_pvalue, DP_pvalue, CSS_test_result, file = file.path(rep_dir, "results.RData"))
  
  cat(paste("Replication", rep_id, "completed and saved in", rep_dir, "\n"))
}



args=commandArgs(trailingOnly = TRUE)
ind_model=as.numeric(args[1])
rep_id = as.numeric(args[2])

model=c(
  'SmallWorld',
  'ScaleFree',
  'Lattice',
  'TreeStar'
)[ind_model]


n=50
p=120

dir_prefix=paste0('Results/')
CSS_param=list(M=100,L=3)


run_experiment(rep_id, n, p, model, dir_prefix, List_CSS_stat = NULL, CSS_param)

  