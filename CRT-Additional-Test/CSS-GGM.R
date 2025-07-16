# GGM co-sufficiency sampling procedures ---------------------------------------
# DO NOT USE LMFIT!

# Generates a uniform vector orthogonal to the given matrix A.
# Arguments:
#   - A: The input matrix
# Returns:
#   - A uniform vector orthogonal to A.
ortho.vec = function(A) {
  # Perform QR decomposition to obtain orthonormal columns
  qr.dec = qr(cbind(A, rnorm(nrow(A))))
  # Generate the orthogonal vector
  return(qr.Q(qr.dec)[, 1+ncol(A)] * sample(c(-1, 1), 1))
}


# Rotates one coordinate.
# Arguments:
#   - x: The data matrix
#   - i: The index of the coordinate to rotate
#   - graph: Adjacent matrix of the graph 
# Returns:
#   - The rotated coordinate.
rotate = function(x, i, graph) {
  # Find neighbors in the graph for the i-th coordinate
  neighbors = setdiff(which(graph[i, ] != 0), i)
  n = nrow(x)
  vec.1 = matrix(1, nr=n, nc=1)
  new.xi = x[, i]
  
  sample.mean = mean(x[, i])
  
  if(length(neighbors) == 0) {
    # If no neighbors, rotate based on the sample mean
    new.xi = sample.mean + ortho.vec(vec.1) * sqrt(sum((x[, i] - sample.mean)^2))
  } else {
    # Otherwise, rotate based on the residuals of the linear regression
    
    # lr.temp = lm(x[, i]~ x[, neighbors] )
    # lr.noise = lm(rnorm(n) ~ x[, neighbors] )
    # new.residual=lr.noise$residuals
    # new.xi = lr.temp$fitted.values +  sqrt(sum(lr.temp$residuals^2)) * new.residual / sqrt(sum(new.residual^2))
    
    
    lr.temp = lm(x[, i]~ x[, neighbors] )
    new.xi = lr.temp$fitted.values +  sqrt(sum(lr.temp$residuals^2)) * ortho.vec(cbind(vec.1, x[, neighbors])) 
    
    
    # lr.temp = lmfit(x[, neighbors], x[, i] - sample.mean ) 
    # new.xi =  (x[, i] - lr.temp$residuals) + ortho.vec(cbind(vec.1, x[, neighbors])) * sqrt(sum(lr.temp$residuals^2))
    
  }
  
  return(new.xi)
}

# Runs one round of the algorithm
# Arguments:
#   - x: The data matrix
#   - orders: The orders in which the coordinates should be rotated
#   - graph: The graph that describes the relationships between variables
# Returns:
#   - The new data matrix after one round of rotation
run = function(x, orders, graph) {
  newx = x
  for(i in orders) {
    newx[, i] = rotate(newx, i, graph)
  }
  return(newx)
}



# Sampling exchangeable copies for a graph
# Arguments:
#   - X: n x p data matrix
#   - graph: Adjacent matrix of the graph G or neighborhood Ni of each i in [p]
#   - M: Number of copies
#   - L: Number of iterations of sampling each chain
#   - I: Permutation of [p] (or a subset of [p])
# Returns:
#   - A list containing M copies of transformed matrices
exchangeable_sampling <- function(X, graph, M=100, L=1, 
                                  I=seq(1:ncol(X)), bReduced = F) {
  
  # Step 1: Start from X and run the previously defined 'run' function 
  # according to the order of I for L times to generate X_hub.
  X_hub <- X
  for (l in 1:L) {
    X_hub <- run(X_hub, I, graph)
  }
  
  # Initialize an empty list to store the M copies
  X_copies <- vector("list", M)
  
  # Step 2: For each copy, start from X_hub and run the 'run' function
  # according to the reversed order of I for L times to generate each X_tilde.
  for (m in 1:M) {
    X_tilde <- X_hub
    for (l in 1:L) {
      X_tilde <- run(X_tilde, rev(I), graph)
    }
    if(bReduced){ 
      X_copies[[m]] <- X_tilde[,I,drop=F]  # no need to store others for CRT
    }else{
      X_copies[[m]] <- X_tilde
    }
  }
  
  # Output: X_copies containing the M transformed matrices
  return(list(X_hub=X_hub, X_copies=X_copies))
}



