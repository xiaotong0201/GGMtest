
# setting -----------------------------------------------------------------

## Band model
# K times s should be smaller than 2
graph_band=function(p,k,s){
  if(k*abs(s)>=2){
    print('Model Error (Band)')
    return(NULL)
  }
  inv.sigma <- diag(p)
  for (i in 1:p){
    for (j in 1:p){
      if(abs(i - j) >= 1 && abs(i - j) <= k){
        inv.sigma[i,j]=s
        inv.sigma[j,i]=s
      }
    }
  }
  return(list(sigma = solve(inv.sigma), inv.sigma = inv.sigma))
}


## Hub Model
## delete an edge with probability q
graph_hub=function(p,hubsize,q=0,noise=0.5){
  
  Omega1 <- matrix(0, nrow = p, ncol = p)
  for (k in 1:(p/hubsize)) {
    i <- hubsize * (k - 1) + 1
    for( j in (i+1):(i+hubsize-1) ){
      if(rbinom(1,1,q)==0){  
        Omega1[i,j]=1
        Omega1[j,i]=1
      }
    }
  }
  diag(Omega1)=colSums(abs(Omega1)) + noise
  
  # lambda_min <- min(eigen(Omega1, only.values = TRUE)$values)
  # diag(Omega1) <- diag(Omega1) + abs(lambda_min) + 0.05
  
  return(list(sigma = solve(Omega1), inv.sigma = Omega1))
}



## Erdos–Renyi model
## Edge connected with probability q
graph_ER=function(p,s,q,epo){
  set.seed(epo)
  inv.sigma <- matrix(0,p,p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      signal=runif(1, min = s/2, max = 3*s/2)
      if (runif(1) < q){
        inv.sigma[i,j]= signal
        inv.sigma[j,i]= signal
      }
    }
  }
  # Make it positive definite as per LaTeX description
  lambda <- min(eigen(inv.sigma)$values)
  inv.sigma = inv.sigma + (abs(lambda) + 0.05) * diag(p)
  return(list(sigma = solve(inv.sigma), inv.sigma = inv.sigma))
}

## marginally, the edge connect with probability q0
graph_ER_null=function(p,q,q0,er_graph){
  inv.sigma <- er_graph$inv.sigma
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if (runif(1) > q0/q){  # delete w.p. 1- q0/q
        inv.sigma[j,i]=inv.sigma[i,j] = 0
      }
    }
  }
  return(list(inv.sigma=inv.sigma))
}

## Helper function 
#  convert a precision matrix to an adjacency matrix for a graph
prec_to_adj <- function(prec_mat,  threshold = 1e-7) {
  adj_mat= 1*(abs(prec_mat) > threshold)
  diag(adj_mat)=0
  
  return(adj_mat)
}


## For reproducing the example in N. Verzelen and F. Villers 2009
graph_Verzelen= function(p, eta) {
  
  
  num_edges <- floor(eta * p * (p - 1) / 2) # tricky! 2023-11-01
  adjacency_matrix <- matrix(0, nrow = p, ncol = p)
  
  edges <- sample(which(lower.tri(adjacency_matrix,diag = F)), num_edges)
  
  adjacency_matrix[edges] <- runif(num_edges, -1, 1)
  
  adjacency_matrix=adjacency_matrix+t(adjacency_matrix)
  
    Pi = generate_partial_corr_matrix(adjacency_matrix) 
    inv.sigma = -Pi; diag(inv.sigma)=1 
    
    ## standardization
    Sigma=solve(inv.sigma)
    diag_s=diag(Sigma)
    Sigma=Sigma/ sqrt(diag_s %o% diag_s)
    ##
    
    
    return(list(sigma = Sigma, inv.sigma = inv.sigma) )
  }

generate_partial_corr_matrix= function(U) {
  diag_values <- apply(abs(U), 2, sum) + 1e-4
  diag(U) <- diag_values
  Pi <- U / sqrt(diag_values %o% diag_values)
  return(Pi)
}

graph_Verzelen_null=function(p,q,eta,prec){
  g_s= 1*(prec!=0)
  
  num_edges_to_delete <- ceiling(q * p * (p - 1) / 2 * eta)  # tricky! 2023-11-01
  
  gs_ind=which( (g_s != 0)&lower.tri(g_s,diag=F), arr.ind = TRUE)
  edges_to_delete <- sample(nrow(gs_ind), min(nrow(gs_ind), num_edges_to_delete))
  g_minus_q_s <- g_s
  g_minus_q_s[gs_ind[edges_to_delete,]] <- 0
  g_minus_q_s = g_minus_q_s* (t(g_minus_q_s!=0))
  return(list(inv.sigma=g_minus_q_s))
  
}


