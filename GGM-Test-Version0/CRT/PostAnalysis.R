
# Load necessary libraries
library(reshape)
library(ggplot2)
source('PostSelection.R')

# Set some variables
bSelectMethod = TRUE  # Set to TRUE if you want to select methods, otherwise FALSE
ind.experiment = 5

# Loop over experiments
for (ind.experiment in 1:8) {
  
  if(ind.experiment==4|ind.experiment==8){
    Max.epo=800
  }else{
    Max.epo=400
  }
  
  # Set the path to the experiment data
  path_prefix = paste0('../../../Programming/LocalExperiments/GGM/CRT/Experiments-', ind.experiment, '/')
  all_output = list()
  
  # Load output data for each epoch
  for (epo in 1:Max.epo) {
    load(file = paste0(path_prefix, epo, '.RData'))
    all_output[[epo]] = output
  }
  
  # check
  sapply(all_output,length)
  
  # Print some information
  print(popul_param$setting)
  
  # Initialize lists to store combined p-values and use time
  combined_pvalues = list()
  combined_usetime = list()
  
  # Loop over signal indices
  for (ind.signal in 1:length(popul_param$modelParam$grid_signal)) {
    combined_pvalues[[ind.signal]] = t(sapply(all_output, function(v) v[[ind.signal]]$pvalues))
    combined_usetime[[ind.signal]] = t(sapply(all_output, function(v) v[[ind.signal]]$use.time))
  }
  
  # check
  sapply(combined_pvalues,dim)
  
  # Look at the time
  matplot(apply(simplify2array(combined_usetime), 3:2, mean))
  
  # Define a function for standard error calculation
  Fun_SE = function(x) {
    x=x[!is.na(x)]
    sd(x) / sqrt(length(x))
  }
  
  # Set significance level
  pvalues_array0 = simplify2array(combined_pvalues)
  sig.lv = 0.05
  
  # Initialize a list to store selected methods
  selected_methods = list()
  
  if (bSelectMethod) {
    # If method selection is enabled, call MethodSelection function
    selected_methods = MethodSelection(popul_param$setting)
    pvalues_array = pvalues_array0[, selected_methods$ind, ]
    print(dimnames(pvalues_array)[[2]])
    print(selected_methods$names)  # check if the two match
    # dimnames(pvalues_array)[[2]] =  selected_methods$names
  } else {
    pvalues_array = pvalues_array0
    selected_methods$names = dimnames(pvalues_array)[[2]]
    selected_methods$theta_lim=c(0, max(popul_param$modelParam$grid_signal))
  }
  
  # Calculate power and standard error tables
  power_table = data.frame(apply(pvalues_array < sig.lv, 3:2, mean,na.rm=T))
  power_SE_table = data.frame(apply(pvalues_array < sig.lv, 3:2, Fun_SE))
  
  
  apply(pvalues_array, 3:2, function(x)sum(!is.na(x)))  # CGM has one NA for theta=3.5
  
  
  print(cbind(popul_param$modelParam$grid_signal, power_table))
  
  num_methods = ncol(power_table)
  power_table$theta = popul_param$modelParam$grid_signal
  power_SE_table$theta = popul_param$modelParam$grid_signal
  
  # Melt the data frames
  df_melt = melt(power_table, id.vars = "theta", variable_name = "Test")
  df_melt$Power = df_melt$value
  df_melt$value = NULL
  
  SE_melt <- melt(power_SE_table, id.vars = "theta", variable_name = "Test")
  SE_melt$SE = SE_melt$value
  SE_melt$value = NULL
  
  # Merge the data and SE data frames
  merged_df <- merge(df_melt, SE_melt, by = c("theta", "Test"), sort = FALSE)
  
  A = dimnames(pvalues_array)[[2]]
  B = selected_methods$names
  merged_df$Test = B[match(merged_df$Test, A)]
  merged_df$Test = factor(merged_df$Test, B)
  
  # Define some aesthetics for ggplot
  distinct_colors <- scales::hue_pal()(num_methods)
  distinct_shapes <- rep(c(1, 2, 3, 4, 5, 6, 8, 15, 16, 17, 18), length.out = num_methods)
  distinct_linewide <- rep(c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6), length.out = num_methods)
  distinct_linetypes <- rep(c("dotdash", "dashed", "solid", "twodash", "longdash", "dotted"), length.out = num_methods)
  
  
  
  # Create ggplot visualization
  gp = ggplot(merged_df, aes(x = theta, y = Power, color = Test, shape = Test, 
                             linetype = Test, linewidth = Test)) +
    geom_line() +
    labs(selected_methods$names) +
    scale_linetype_manual(values = distinct_linetypes) +
    scale_linewidth_manual(values = distinct_linewide) +
    scale_color_manual(values = distinct_colors) +
    theme_bw() +
    theme(legend.position = "right",
          legend.key.size = unit(1.2, 'cm'), 
          legend.key.width = unit(1.5, 'cm'),
          text = element_text(size = 20), # Adjusts overall text size; you can change the value as needed
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.title.x = element_text(size = 20), # Adjusts x-axis label size
          axis.title.y = element_text(size = 20), # Adjusts y-axis label size
          plot.title = element_text(size = 20) # Adjusts plot title size
    ) +
    geom_point(size = 3) +
    labs(x = expression(theta),
         y = "Power")  +
    xlim(selected_methods$theta_lim[1] - 0.1, selected_methods$theta_lim[2] + 0.1) + 
    ylim(0,     min(   max(merged_df$Power+merged_df$SE) * ifelse(ind.experiment==4,1, 1.2) , 1)  ) + 
    scale_shape_manual(values = distinct_shapes) +
    scale_fill_manual(values = distinct_colors) +
    geom_errorbar(aes(ymin = Power - SE, ymax = Power + SE), 
                  linetype = "solid", width = 0.1, linewidth = 0.4,
                  show.legend = FALSE)
  
  print(gp)
  
  
  new_gp <-  gp +
    geom_hline(aes(yintercept = 0.05, size = "Level = .05"), color = "black", linetype = "dotted") +
    scale_size_manual(name = "Parameters", values = c(0.5, "Level = .05" = 0.5)) +
    guides(size = guide_legend(title = "\n "))
  
  
  # Print the plot
  print(new_gp)
  
  
  
  # Save the plot as an EPS file
  ggsave(filename = paste0('CRT-', popul_param$setting, ".pdf"), plot = new_gp, device = "pdf", width = 8, height = 6)
}
