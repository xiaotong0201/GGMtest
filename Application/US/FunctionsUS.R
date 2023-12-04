
# Loading required libraries
library(Matrix)
library(MASS)
library(mvtnorm)
library(glasso)
library(Rfast)
library(glmnet)
library(rags2ridges)
source('../../CSS-GGM.R') # Source external scripts
source('../../GGMTesting.R')
source('../../GoF/GoF-Statistics.R')

# Function to run all goodness-of-fit tests
run_all_gof_tests <- function(x, G0, List_CSS_stat = NULL, CSS_param = list(M = 100, L = 3),
                              CSS = NULL, dir_prefix = '', bSave=F, bReturnCSS=F) {
  if(bSave){
    # Path setup for data storage
    path_prefix <- dir_prefix
    if (!file.exists(path_prefix)) {
      dir.create(path_prefix) # Create directory if it doesn't exist
    }
    file_prefix <- paste0('US_n=', nrow(x))
    file_prefix <- paste0(path_prefix, file_prefix)
  }
  
  # Goodness-of-fit tests
  VV_pvalue <- VV_GoF(x, G0) # Verzelen et al
  DP_pvalue <- Bonf_GoF(x, G0) # Bonferroni
  
  print( paste0('VV:',VV_pvalue))
  print( paste0('DP:',DP_pvalue))
  
  # CSS sampling
  if (is.null(CSS)) {
    start.time <- Sys.time()
    cat('CSS begin;')
    CSS <- exchangeable_sampling(x, G0, M = CSS_param$M, L = CSS_param$L)
    time.taken <- Sys.time() - start.time
    cat('Time used: ', time.taken, 's;')
    cat('CSS end;\n')
    if(bSave){
      save(CSS, file = paste0(file_prefix, "_CSS", ".RData"))
    }
  }
  
  # Calculate statistics
  cat('CSS testing begin;')
  CSS_test_result <- list()
  for (i in seq_along(List_CSS_stat)) {
    f <- List_CSS_stat[i]
    cat(f);cat(':')
    start.time <- Sys.time()
    CSS_test_result[[i]] <- CSSGoF(x, G0, CSS$X_copies, Fun_Stat = get(f))
    print(CSS_test_result[[i]]$pvalue)
    time.taken <- Sys.time() - start.time
    cat('Time used: ', time.taken, 's;')
    CSS_test_result[[i]]$time.taken <- time.taken
  }
  cat('CSS testing end;\n')
  
  if(bSave){
    save(CSS_test_result, file = paste0(file_prefix, "_TestStat", ".RData"))
  }
  
  CSS_test_pvalues <- sapply(CSS_test_result, `[[`, "pvalue")
  names(CSS_test_pvalues) <- List_CSS_stat
  
  # Wrapping results in a list
  final_result <- list(
    VV_pvalue = VV_pvalue,
    DP_pvalue = DP_pvalue,
    CSS_test_pvalues = CSS_test_pvalues
  )
  
  if(bReturnCSS){
    final_result$CSS=CSS
  }
  return(final_result)
}
