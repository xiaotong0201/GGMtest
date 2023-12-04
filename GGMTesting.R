#-----------------------------------------------------------
# Function: ComputePValue
# Purpose:  Compute the p-value based on observed and sampled test statistics
# Arguments:
#   T0       - Observed test statistic
#   Ttilde   - Sampled test statistics
#   type     - Type of test ('One-sided' or 'Two-sided')
#   randomize- Whether to randomize in case of ties
# Returns:   Computed p-value
#-----------------------------------------------------------

ComputePValue <- function(T0, Ttilde, type =  "One-sided", #c('One-sided', 'Two-sided')
                          randomize = TRUE) {
  
  
  if(is.matrix(Ttilde)){
    cntGreater <- rowSums(Ttilde > T0)
    cntSmaller <- rowSums(Ttilde < T0)
    cntEqual   <- rowSums(Ttilde == T0) + 1
    M <-  ncol(Ttilde)
  }else{
    
    cntGreater <- sum(Ttilde > T0)
    cntSmaller <- sum(Ttilde < T0)
    cntEqual   <- sum(Ttilde == T0) + 1
     M <- length(Ttilde)
  }
  
  # Random sampling for tie-breaking
  S <- ifelse(randomize, 
              sapply(cntEqual,function(x)sample(x,1)) # cntEqual may be a vector
              , cntEqual)
  
  # Calculate p-value
  if (type == "One-sided") {
    pvalue <- (cntGreater + S) / (1 + M)
  } else {
    pvalue <- 2 * ( pmin(cntGreater, cntSmaller) + S) / (1 + M)
  }
  
  return(pvalue)
}

#-----------------------------------------------------------
# Function: CSSGoF
# Purpose:  Conducts a goodness-of-fit test based on exchangeable sampling
# Arguments:
#   X        - Data matrix
#   graph    - Adjacency matrix of the graph
#   X_copies - Optional, precomputed exchangeable samples
#   Fun_Stat - Function to compute test statistic
#   ...      - Additional arguments for exchangeable_sampling function
# Returns:   Computed p-value
#-----------------------------------------------------------
CSSGoF <- function(X, graph, X_copies = NULL, 
                   Fun_Stat, type='One-sided', randomize=T,M=100,L=1) {
  
  # Generate exchangeable samples if not provided
  if (is.null(X_copies)) {
    CSS_samples <- exchangeable_sampling(X, graph, M=M, L=L)
    X_copies <- CSS_samples$X_copies
  }
  
  T0 <- Fun_Stat(X, graph)
  Ttilde <- simplify2array(lapply(X_copies, Fun_Stat, graph))
  
  pvalue=ComputePValue(T0, Ttilde, type=type, randomize = randomize)
  return(list(pvalue=pvalue, T0=T0, Ttilde=Ttilde))
}
