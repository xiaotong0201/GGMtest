
CRT_oversized_graphs = function(popul_param,
               CSS_param = list(M = 100, L = 3),
               epo, 
               CSS_sample = NULL, dir_prefix = '') {
  
  # Extract parameters
  p <- popul_param$P
  n <- popul_param$N
  setT <- popul_param$setT
  setting <- popul_param$setting
  modelparam <- popul_param$modelParam
  
  # Print initial parameters
  print(paste(n, p, setting, sep = ","))
  
  # Create directory for data storage
  path_prefix = paste0(dir_prefix, 'DataStorage-', format(Sys.time(), "%b-%e-%Y"), '/')
  if (!file.exists(path_prefix)) {
    dir.create(path_prefix)
  }
  
  # Construct file prefix
  file_prefix = paste0(path_prefix, paste0(c(rbind(names(popul_param)[1:3], unlist(popul_param[1:3]))), collapse = "_"))
  
  # Generate graphs
  true_graph = graph_band(p, modelparam$K, modelparam$s)
  
  
  gap_bandwidth=ifelse(setting=='l_l_w', 2, 6)
  oversize_graph_1 = graph_band(p, modelparam$K + gap_bandwidth*1, 1e-3)
  oversize_graph_2 = graph_band(p, modelparam$K + gap_bandwidth*2, 1e-3)
  
  # Permute variables and calculate sigma
  set.seed(epo)
  perm.var = sample(p, p)
  sigma = true_graph$sigma[perm.var, perm.var]
  inv.sigma = true_graph$inv.sigma[perm.var, perm.var]
  
  # Normalize sigma and inv.sigma
  diag_s = sqrt(diag(sigma))
  sigma = sweep(sweep(sigma, 1, 1 / diag_s, "*"), 2, 1 / diag_s, "*")
  inv.sigma = sweep(sweep(inv.sigma, 1, diag_s, "*"), 2, diag_s, "*")
  
  # Convert precision matrix to adjacency matrix
  G = prec_to_adj(inv.sigma)
  
  # Check for model errors
  if (min(eigen(inv.sigma)$values) < 0) {
    print('Model Error!')
    return(NULL)
  }
  if (2 * max(colSums(G)) + 3 > n) {
    print('True graph too large!')
    return(NULL)
  }
  
  # List of graphs for CSS
  Graphs_CSS = list(G, prec_to_adj(oversize_graph_1$inv.sigma[perm.var, perm.var]), prec_to_adj(oversize_graph_2$inv.sigma[perm.var, perm.var]))
  
  # Simulate data
  x = simulate_data(sigma, n)
  beta = random_beta(p)
  noise_unif = runif(n)
  
  # CSS sampling and testing
  results_graphsize = list()
  for (i_graph in seq_along(Graphs_CSS)) {
    start.time = Sys.time()
    cat('CSS begin;')
    CSS_sample = exchangeable_sampling(x, Graphs_CSS[[i_graph]], I = setT, bReduced = TRUE, M = CSS_param$M, L = CSS_param$L)
    time.taken = Sys.time() - start.time
    cat('Time used: ', time.taken, 's;')
    cat('CSS end;\n')
    
    List_Stat_Comb = GetStatComb(setting)
    test_results = list()
    
    for (ind_theta in seq_along(modelparam$grid_signal)) {
      theta = modelparam$grid_signal[ind_theta]
      print(theta)
      
      y = simulate_y(x, setting, setT, theta, beta, noise_unif = noise_unif)
      
      # Run tests and store results
      new_test_result = Run_all_test(x, y, setT, CSS_sample, setting, List_Stat_Comb,IncludeBaseline = F)
      test_results[[ind_theta]] = list(
        theta = theta,
        pvalues = new_test_result$pvalues,
        use.time = new_test_result$use.time
      )
      print(unlist(test_results[[ind_theta]]))
    }
    results_graphsize[[i_graph]] = test_results
  }
  
  return(results_graphsize)
}




