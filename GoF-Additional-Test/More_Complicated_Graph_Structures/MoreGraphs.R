library(igraph)

# Function to create a ring lattice
create_ring_lattice <- function(n, k) {
  if (k %% 2 != 0) stop("k must be an even number.")
  
  adj_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:(k / 2)) {
      adj_matrix[i, ((i - 1 + j) %% n) + 1] <- 1
      adj_matrix[((i - 1 + j) %% n) + 1, i] <- 1
    }
  }
  return(adj_matrix)
}

# Function to perform rewiring on an adjacency matrix
rewire_graph <- function(adj_matrix, prob) {
  n <- nrow(adj_matrix)
  k <- sum(adj_matrix[1, ]) # Degree of each node (assuming uniform degree)
  
  for (distance in 1:(k / 2)) {
    for (i in 1:n) {
      if (runif(1) < prob) {
        # Find a new target node to rewire
        possible_targets <- setdiff(1:n, c(i, which(adj_matrix[i, ] == 1)))
        if (length(possible_targets) > 0) {
          new_target <- sample(possible_targets, 1)
          
          # Rewire edge
          old_target <- ((i - 1 + distance) %% n) + 1
          adj_matrix[i, old_target] <- 0
          adj_matrix[old_target, i] <- 0
          
          adj_matrix[i, new_target] <- 1
          adj_matrix[new_target, i] <- 1
        }
      }
    }
  }
  
  return(adj_matrix)
}

# Function to generate a scale-free graph
generate_scale_free_graph <- function(p, m, m0) {
  if (m > m0) stop("m0 must be greater than or equal to m")
  if (p < m0) stop("p must be greater than or equal to m0")
  
  adj_matrix <- matrix(0, nrow = p, ncol = p)
  adj_matrix[1:m0, 1:m0] <- 1
  diag(adj_matrix) <- 0 # No self-loops
  
  for (new_node in (m0 + 1):p) {
    degrees <- colSums(adj_matrix[1:(new_node - 1), 1:(new_node - 1)])
    prob <- degrees / sum(degrees)
    connections <- sample(1:(new_node - 1), size = m, prob = prob, replace = FALSE)
    adj_matrix[new_node, connections] <- 1
    adj_matrix[connections, new_node] <- 1
  }
  
  return(adj_matrix)
}

# Function to simulate node addition/removal
update_graph <- function(G, q_del, q_add) {
  p <- nrow(G)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      if (G[i, j] == 1) {
        if (runif(1) < q_del) G[j, i] <- G[i, j] <- 0
      } else {
        if (runif(1) < q_add) G[j, i] <- G[i, j] <- 1
      }
    }
  }
  return(G)
}

# Generate random signals for weights
generate_random_signals <- function(p) {
  random_signals <- matrix(sample(c(-1, 1), p * p, replace = TRUE), p)
  random_signals[upper.tri(random_signals, diag = TRUE)] <- 0
  random_signals <- random_signals + t(random_signals)
  return(random_signals)
}

# Adjust Omega to positive definite Sigma
adjust_omega_to_sigma <- function(Omega, constant_base) {
  diag(Omega) <- constant_base + abs(min(eigen(Omega)$values))
  Sigma <- cov2cor(solve(Omega))
  if (min(eigen(Sigma)$values) <= 0) {
    stop("Generated Sigma is not positive definite.")
  }
  return(Sigma)
}

# Setup graph and generate Sigma for a given model
setup_graph_and_sigma <- function(model, p, q_del = 0.2, q_add = 0.1, 
                                  constant_base = 2) {
  if (model == "SmallWorld") {
    G0 <- create_ring_lattice(p,8)
    G0 <- rewire_graph(G0, 0.2)
    graph_true <- update_graph(G0, q_del, q_add)
    random_signals <- generate_random_signals(p)
    Omega <- graph_true * random_signals
    Sigma <- adjust_omega_to_sigma(Omega, constant_base)
  } else if (model == "ScaleFree") {
    G0 <- generate_scale_free_graph(p, m = 2, m0 = 5)
    graph_true <- update_graph(G0, q_del, q_add)
    random_signals <- generate_random_signals(p)
    Omega <- graph_true * random_signals
    Sigma <- adjust_omega_to_sigma(Omega, constant_base )
  } else if (model == "Lattice") {
    G0 <- as.matrix(get.adjacency(make_lattice(c(floor(sqrt(p)), ceiling(sqrt(p)) + 1), nei = 1)))
    graph_true <- as.matrix(get.adjacency(make_lattice(c(floor(sqrt(p)), ceiling(sqrt(p)) + 1), nei = 2)))
    random_signals <- generate_random_signals(p)
    Omega <- graph_true * random_signals
    Sigma <- adjust_omega_to_sigma(Omega, constant_base )
  } else if (model == "StarTree") {
    G0 <- matrix(0, p, p)
    G0[1, 2:p] <- 1
    G0[2:p, 1] <- 1
    graph_true <- as.matrix(get.adjacency(make_tree(p, children = 3, mode = "undirected")))
    random_signals <- generate_random_signals(p)
    Omega <- graph_true * random_signals
    Sigma <- adjust_omega_to_sigma(Omega, constant_base )
  } else if(model=='TreeStar'){
    G0 <- as.matrix(get.adjacency(make_tree(p, children = 3, mode = "undirected")))
    graph_true <- matrix(0, p, p)
    graph_true[1, 2:p] <- 1
    graph_true[2:p, 1] <- 1
    
    random_signals <- generate_random_signals(p)
    Omega <- graph_true * random_signals
    Sigma <- adjust_omega_to_sigma(Omega, constant_base )
  } else {
    stop("Invalid model specified.")
  }
  
  return(list(G0 = G0, Sigma = Sigma))
}
