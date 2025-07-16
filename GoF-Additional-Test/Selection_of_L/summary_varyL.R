## Load necessary libraries
library(dplyr)
library(ggplot2)

## Parameters
param_DenseER_template <- list(
  model = "ER",
  N = c(50),
  P = c(120),
  modelParam = list(
    Q = c(0.5),
    Q0 = 0.25,
    s = c(0.01, 0.02)
  )
)

List_CSS_stat = c(
  "PRC_SS", "PRC_SA", "PRC_Z_SS", "PRC_Z_SA",
  "ERC_SS", "ERC_SA", "ERC_Z_SS", "ERC_Z_SA",
  "F_max", "F_sum", "F_logsum", "F_Z_SA", "F_Z_SS", 
  "glr_glasso"
)

L_values <- c(1, 2, 3, 5, 10, 20)  # Values of L to test
sig.lv <- 0.05  # Significance level

read_experiments_epo <- function(param_template, output_list, L_values, List_CSS_stat, sig.lv = 0.05) {
  MaxEpo <- length(output_list)
  summary_result <- list()
  lFR <- 0
  
  for (ind.N in 1:length(param_template$N)) {
    N <- param_template$N[ind.N]
    for (ind.P in 1:length(param_template$P)) {
      P <- param_template$P[ind.P]
      
      combos <- expand.grid(param_template$modelParam)
      
      # Loop over each combination of model parameters
      for (ind.setup in 1:nrow(combos)) {
        np.id <- (ind.N - 1) * length(param_template$P) + ind.P
        exp.id <- ind.setup + nrow(combos) * (np.id - 1)
        
        combined_results <- list()
        VV_pvalues <- c()  # Store VV p-values for each epoch
        DP_pvalues <- c()  # Store DP p-values for each epoch
        
        
        for (epo in 1:MaxEpo) {
          css_results <- output_list[[epo]][[exp.id]]$CSS_test_results
          VV_pvalues <- c(VV_pvalues, output_list[[epo]][[exp.id]]$VV_pvalue)
          DP_pvalues <- c(DP_pvalues, output_list[[epo]][[exp.id]]$DP_pvalue)
          
          
          for (L in L_values) {
            L_key <- paste0("L=", L)
            if (!is.null(css_results[[L_key]])) {
              stat_pvalues <- c()  # Use a named vector to store p-values
              
              for (stat_key in 1:length(List_CSS_stat)) {
                stat_name = List_CSS_stat[stat_key]
                if (!is.null(css_results[[L_key]][[stat_key]]$pvalue)) {
                  stat_pvalues[stat_name] <- css_results[[L_key]][[stat_key]]$pvalue
                } else {
                  stat_pvalues[stat_name] <- NA  # Fill missing values with NA
                }
              }
              
              # Combine all statistics for this `L` into a single row
              combined_results <- bind_rows(
                combined_results,
                cbind(data.frame(L = L, epo = epo), as.data.frame(t(stat_pvalues)))
              )
            }
          }
          
          
        }
        combined_results <- as.data.frame(combined_results)
        temp=subset(combined_results, L==1);colMeans(temp<0.05)
        temp=subset(combined_results, L==20);colMeans(temp<0.05)
        
        
        temp1=subset(combined_results, L==1)
        temp2=subset(combined_results, L==20)
        plot(temp1[,'F_sum'],temp2[,'F_sum'])
        
        # Summarize results by `L`
        if (nrow(combined_results) > 0) {
          # Compute power and SE for each statistic
          power <- combined_results %>%
            group_by(L) %>%
            summarize(across(where(is.numeric), ~ mean(. < sig.lv, na.rm = TRUE), .names = "power_{col}"))
          
          SE <- combined_results %>%
            group_by(L) %>%
            summarize(across(where(is.numeric), ~ sd(. < sig.lv, na.rm = TRUE) / sqrt(n()), .names = "SE_{col}"))
        } else {
          power <- data.frame(L = L_values, statistic = NA, power = NA)
          SE <- data.frame(L = L_values, statistic = NA, SE = NA)
        }
        
        
        # Summarize VV and DP p-values across epochs
        VV_summary <- data.frame(
          statistic = "VV_pvalue",
          mean_pvalue = mean(VV_pvalues<sig.lv, na.rm = TRUE),
          SE_pvalue = sd(VV_pvalues<sig.lv, na.rm = TRUE) / sqrt(length(VV_pvalues))
        )
        
        DP_summary <- data.frame(
          statistic = "DP_pvalue",
          mean_pvalue = mean(DP_pvalues<sig.lv, na.rm = TRUE),
          SE_pvalue = sd(DP_pvalues<sig.lv, na.rm = TRUE) / sqrt(length(DP_pvalues))
        )
        
        modelParam <- as.list(combos[ind.setup, ])
        names(modelParam) <- names(param_template$modelParam)
        
        summary_result[[lFR + 1]] <- list(
          model = param_template$model,
          N = N,
          P = P,
          modelParam = modelParam,
          combined_results = combined_results,
          power = power,
          SE = SE,
          VV_summary = VV_summary,
          DP_summary = DP_summary
        )
        lFR <- lFR + 1
      }
    }
  }
  return(summary_result)
}


## Function to summarize experiments
summarize_experiments <- function(base_path, epo_range, L_values) {
  output_list <- list()
  
  # Loop through each epo and load the corresponding .RData file
  for (epo in epo_range) {
    file_path <- paste0(base_path, "Experiments-5/", epo, ".RData")
    if (file.exists(file_path)) {
      load(file_path)  # Assumes `output` is the object in the .RData file
      output_list[[epo]] <- output
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  
  # Read and summarize the experiments
  summary_result <- read_experiments_epo(param_DenseER_template, output_list, L_values, List_CSS_stat, sig.lv)
  return(summary_result)
}

## Base path and settings
base_path <- "../../../Programming/LocalExperiments/GGM/GoF/"
epo_range <- 1:400  # Range of epo to process

## Summarize results
summary_result <- summarize_experiments(base_path, epo_range, L_values)

View(summary_result[[1]]$power)
summary_result[[1]]$VV_summary
summary_result[[1]]$DP_summary

View(summary_result[[2]]$power)
summary_result[[2]]$VV_summary
summary_result[[2]]$DP_summary


## Example: Extract power metric
power_table <- extract_metric("power", summary_result)
