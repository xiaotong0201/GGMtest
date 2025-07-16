load_all_data <- function(results_dir = "Results", sig.lv = 0.05, List_CSS_stat = c(
  "PRC_SS", "PRC_SA", "PRC_Z_SS", "PRC_Z_SA",
  "ERC_SS", "ERC_SA", "ERC_Z_SS", "ERC_Z_SA",
  "F_max", "F_sum", "glr_glasso"
)) {
  # List all RData files recursively
  rdata_files <- list.files(results_dir, pattern = "results.RData", recursive = TRUE, full.names = TRUE)
  
  if (length(rdata_files) == 0) {
    stop("No results.RData files found in the specified directory.")
  }
  
  # Initialize storage
  output_list <- list()
  
  # Load all data
  for (file in rdata_files) {
    load(file)  # Loads VV_pvalue, DP_pvalue, CSS_test_result into the environment
    
    # Extract model and replication ID from file path
    path_parts <- strsplit(file, "/")[[1]]
    model <- path_parts[2]
    replication <- as.numeric(gsub("rep_", "", path_parts[3]))
    
    # Store results in output_list
    if (!is.list(output_list[[model]])) {
      output_list[[model]] <- list()
    }
    output_list[[model]][[replication]] <- list(
      VV_pvalue = VV_pvalue,
      DP_pvalue = DP_pvalue,
      CSS_test_pvalues = sapply(CSS_test_result, function(res) res$pvalue)
    )
  }
  
  # Prepare a combined p-value table
  pvalue_tables <- list()
  for (model in names(output_list)) {
    pvalue.table <- do.call(rbind, lapply(output_list[[model]], function(res) {
      c(res$VV_pvalue, res$DP_pvalue, unlist(res$CSS_test_pvalues))
    }))
    colnames(pvalue.table) <- c("VV", "DP", List_CSS_stat)
    pvalue_tables[[model]] <- pvalue.table
  }
  
  return(list(pvalue_tables = pvalue_tables, output_list = output_list))
}

summarize_data <- function(pvalue_tables, sig.lv = 0.05) {
  summary_result <- list()
  lFR <- 0
  
  for (model in names(pvalue_tables)) {
    pvalue.table <- pvalue_tables[[model]]
    
    # Compute power and standard error
    power <- colMeans(pvalue.table < sig.lv)
    SE <- apply(pvalue.table < sig.lv, 2, function(v) sqrt(var(v) / length(v)))
    
    # Summarize results
    summary_result[[lFR + 1]] <- list(
      model = model,
      pvalue.table = pvalue.table,
      power = power,
      SE = SE
    )
    lFR <- lFR + 1
  }
  
  return(summary_result)
}
create_power_table <- function(summary_result, List_CSS_stat = c(
  "PRC_SS", "PRC_SA", "PRC_Z_SS", "PRC_Z_SA",
  "ERC_SS", "ERC_SA", "ERC_Z_SS", "ERC_Z_SA",
  "F_max", "F_sum", "glr_glasso"
)) {
  # Include "VV" and "DP" in the list of methods
  methods <- c("VV", "DP", List_CSS_stat)
  
  # Initialize a table with models as rows and methods as columns
  models <- sapply(summary_result, function(res) res$model)
  power_table <- matrix(NA, nrow = length(models), ncol = length(methods),
                        dimnames = list(models, methods))
  
  # Fill the power table with values
  for (i in seq_along(summary_result)) {
    power <- summary_result[[i]]$power
    power_table[i, ] <- power
  }
  
  return(power_table)
}


create_se_table <- function(summary_result, List_CSS_stat = c(
  "PRC_SS", "PRC_SA", "PRC_Z_SS", "PRC_Z_SA",
  "ERC_SS", "ERC_SA", "ERC_Z_SS", "ERC_Z_SA",
  "F_max", "F_sum", "glr_glasso"
)) {
  # Include "VV" and "DP" in the list of methods
  methods <- c("VV", "DP", List_CSS_stat)
  
  # Initialize a table with models as rows and methods as columns
  models <- sapply(summary_result, function(res) res$model)
  se_table <- matrix(NA, nrow = length(models), ncol = length(methods),
                        dimnames = list(models, methods))
  
  # Fill the se table with values
  for (i in seq_along(summary_result)) {
    se <- summary_result[[i]]$SE
    se_table[i, ] <- se
  }
  
  return(se_table)
}

# Step 1: Load all data
loaded_data <- load_all_data("Results")

# Step 2: Summarize the data
summary_result <- summarize_data(loaded_data$pvalue_tables)

# Example: Print the summary for the first model
print(summary_result[[1]]$power)
print(summary_result[[2]]$power)
print(summary_result[[3]]$power)
print(summary_result[[4]]$power)

# Step 3

# Example usage
power_table <- create_power_table(summary_result)

# Print the power table
print(power_table)


se_table <- create_se_table(summary_result)
print(se_table)



###############

library(xtable)

# Function to create a filtered and reordered LaTeX table
create_filtered_and_ordered_latex_table <- function(
    power_table,
    se_table,
    selected_methods = c("VV", "DP", "PRC_SS", "ERC_SS", "F_sum", "glr_glasso"),
    custom_names = c("VV", "DP", "PRC", "ERC", "F$_\\Sigma$", "GLR-$\\ell_1$"),
    row_order = c("StarTree", "Lattice", "SmallWorld", "ScaleFree"),
    caption = "Power Comparison for Selected Methods",
    label = "tab:filtered_power_table",
    output_file = "filtered_power_table.tex"
) {  # Subset the tables for selected methods
  filtered_power_table <- power_table[, selected_methods, drop = FALSE]
  filtered_se_table <- se_table[, selected_methods, drop = FALSE]
  colnames(filtered_power_table) <- custom_names
  colnames(filtered_se_table) <- custom_names
  
  # Reorder rows
  filtered_power_table <- filtered_power_table[row_order, , drop = FALSE]
  filtered_se_table <- filtered_se_table[row_order, , drop = FALSE]
  
  # Combine power and SE into a single table with "Power (SE)" format
  combined_table <- apply(filtered_power_table, c(1, 2), as.character)
  for (i in seq_len(nrow(filtered_power_table))) {
    for (j in seq_len(ncol(filtered_power_table))) {
      combined_table[i, j] <- sprintf("%.3f (%.3f)", filtered_power_table[i, j], filtered_se_table[i, j])
    }
  }
  
  # Convert to a data frame and add model names
  df_combined <- as.data.frame(combined_table, stringsAsFactors = FALSE)
  df_combined <- cbind(Model = rownames(filtered_power_table), df_combined)
  
  # Generate LaTeX code using xtable
  latex_table <- xtable(
    df_combined,
    caption = caption,
    label = label,
    align = c("l", "l", rep("c", ncol(filtered_power_table)))  # Align: left for first column, center for others
  )
  
  # Save LaTeX table to a file
  sink(output_file)
  print(
    latex_table,
    include.rownames = FALSE,  # Do not include row numbers
    sanitize.text.function = identity,  # Keep LaTeX formatting intact
    table.placement = "ht"  # Place the table here or at the top of the page
  )
  sink()
  
  message("LaTeX table saved to ", output_file)
  }

# Example usage
create_filtered_and_ordered_latex_table(
  power_table,  se_table,
  selected_methods = c("VV", "DP", "PRC_SS", "ERC_SS", "F_sum", "glr_glasso"),
  custom_names = c("$M^1P_1$", "Bonf", "PRC", "ERC", "F$_\\Sigma$", "GLR-$\\ell_1$"),
  row_order = c("TreeStar", "Lattice", "SmallWorld", "ScaleFree"),
  caption = "Power Comparison for Selected Methods",
  label = "tab:filtered_power_table",
  output_file = "filtered_power_table.tex"
)



