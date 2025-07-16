### Completed simulation ( 31 Oct 2023) should be found in backup
### This version (2023-11-04) focus on only a few stat

# Load necessary libraries
library(reshape)
library(ggplot2)
source('PostSelection.R')

# Loop over selected experiments
for (ind.experiment in c(1, 5)) {
  
  
  
  if(ind.experiment==4|ind.experiment==8){
    Max.epo=800
  }else{
    Max.epo=400
  }
  
  # Set the path to the experiment data
  path_prefix = paste0('/Users/stahd/GitHub/Programming/LocalExperiments/GGM/CRT-More/Experiments-', ind.experiment, '/')
  all_output = list()
  
  # Load output data for each epoch
  for (epo in 1:Max.epo) {
    load(file = paste0(path_prefix, epo, '.RData'))
    all_output[[epo]] = output
  }
  
  # Check the length of all_output
  sapply(all_output, length)
  
  # Print some information
  print(popul_param$setting)
  
  array_pvalues=list()
  array_usetime=list()
  
  for (ind.graph in 1:3){
    # Initialize lists to store combined p-values and use time
    combined_pvalues = list()
    combined_usetime = list()
    
    # Loop over signal indices
    for (ind.signal in 1:length(popul_param$modelParam$grid_signal)) {
      combined_pvalues[[ind.signal]] = t(sapply(all_output, function(v) v[[ind.graph]][[ind.signal]]$pvalues))
      combined_usetime[[ind.signal]] = t(sapply(all_output, function(v) v[[ind.graph]][[ind.signal]]$use.time))
    }
    
    
    
    # Check dimensions of combined_pvalues
    sapply(combined_pvalues, dim)
    
    # Look at the time
    matplot(apply(simplify2array(combined_usetime), 3:2, mean))
    
    array_pvalues[[ind.graph]]=combined_pvalues
    array_usetime[[ind.graph]]=combined_usetime
  }
  
  ### below is left 
  
  # Define a function for standard error calculation
  Fun_SE = function(x) {
    x = x[!is.na(x)]
    sd(x) / sqrt(length(x))
  }
  
  # Set significance level
  sig.lv = 0.05
  
  
  # Initialize a list to store selected methods
  selected_methods = list()
  
  # Combine data for all graphs
  combined_df_list = list()
  
  for (ind.graph in 1:3) {
    pvalues_array0 = simplify2array(array_pvalues[[ind.graph]])
    
    
    pvalues_array = pvalues_array0
    selected_methods$names = dimnames(pvalues_array)[[2]]
    selected_methods$theta_lim = c(0, max(popul_param$modelParam$grid_signal))
    
    
    
    # Calculate power and standard error tables
    power_table = data.frame(apply(pvalues_array < sig.lv, 3:2, mean, na.rm = TRUE))
    power_SE_table = data.frame(apply(pvalues_array < sig.lv, 3:2, Fun_SE))
    
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
    merged_df$Graph = factor(ind.graph, levels = 1:3, labels = c("Faithful Graph", "Oversized Graph 1", "Oversized Graph 2"))
    
    combined_df_list[[ind.graph]] = merged_df
  }
  
  # Combine all data frames into one
  combined_df = do.call(rbind, combined_df_list)
  
  # Define some aesthetics for ggplot
  distinct_colors <- scales::hue_pal()(length(unique(combined_df$Graph)))
  distinct_shapes <- rep(c(1, 2, 3, 4, 5, 6, 8, 15, 16, 17, 18), length.out = length(unique(combined_df$Graph)))
  distinct_linewide <- rep(c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6), length.out = length(unique(combined_df$Graph)))
  distinct_linetypes <- rep(c("dotdash", "dashed", "solid", "twodash", "longdash", "dotted"), length.out = length(unique(combined_df$Graph)))
  
  # Create ggplot visualization for each test across all graphs
  for (test_name in unique(combined_df$Test)) {
    test_df <- subset(combined_df, Test == test_name)
    
    gp = ggplot(test_df, aes(x = theta, y = Power, color = Graph, shape = Graph, 
                               linetype = Graph, linewidth = Graph, group = Graph)) +
      geom_line() +
      # labs(title = test_name) +
      scale_linetype_manual(values = distinct_linetypes) +
      scale_linewidth_manual(values = distinct_linewide) +
      scale_color_manual(values = distinct_colors) +
      theme_bw() +
      theme(legend.position = "right",
            legend.key.size = unit(1.2, 'cm'), 
            legend.key.width = unit(1.5, 'cm'),
            text = element_text(size = 20),
            axis.text = element_text(size = 20),
            legend.text = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            plot.title = element_text(size = 20)
      ) +
      geom_point(size = 3) +
      labs(x = expression(theta),
           y = "Power")  +
      xlim(selected_methods$theta_lim[1] - 0.1, selected_methods$theta_lim[2] + 0.1) + 
      ylim(0, min(max(test_df$Power + test_df$SE) * ifelse(ind.experiment == 4, 1, 1.2), 1)) + 
      scale_shape_manual(values = distinct_shapes) +
      scale_fill_manual(values = distinct_colors) +
      geom_errorbar(aes(ymin = Power - SE, ymax = Power + SE), 
                    linetype = "solid", width = 0.1, linewidth = 0.4,
                    show.legend = FALSE)
    
    new_gp <- gp +
      geom_hline(aes(yintercept = 0.05, size = "Level = .05"), color = "black", linetype = "dotted") +
      scale_size_manual(name = "Parameters", values = c(0.5, "Level = .05" = 0.5)) +
      guides(size = guide_legend(title = "\n "))
    
    # Print the plot
    print(new_gp)
    
    # Save the plot as a PDF file
    ggsave(filename = paste0('CRT-', popul_param$setting, "-", test_name, "-CombinedGraphs.pdf"), plot = new_gp, device = "pdf", width = 8, height = 6)
  }
}