library(nortest)

# Transformation function
transform_piece_data <- function(x, lambda, cutoff=0.18) {
  c_value <- abs(cutoff)^lambda - 0.01*log(abs(cutoff))
  sapply(x, function(xi) {
    if ( abs(xi) < cutoff) {
      return(sign(xi)*abs(xi)^lambda)
    } else {
      return(sign(xi)*(c_value + 0.01* log(abs(xi))))
    }
  })
}

# Function to test normality and return a p-value
test_normality <- function(data) {
  test <- shapiro.test(data)
  return(test$p.value)
}

# Optimization function to find best lambda and gamma
find_best_parameters <- function(x, range_lambda, range_gamma) {
  best_params <- c(lambda = NA, gamma = NA)
  best_p_value <- 0
  
  lambda_grid <- seq(range_lambda[1], range_lambda[2], length.out = 20)
  gamma_grid <- seq(range_gamma[1], range_gamma[2], length.out = 20)
  
  for (lambda in lambda_grid) {
    for (gamma in gamma_grid) {
      transformed_x <- transform_piece_data(x, lambda, cutoff=gamma)
      p_value <- test_normality(transformed_x)
      
      if (p_value > best_p_value) {
        best_p_value <- p_value
        best_params <- c(lambda, gamma)
      }
    }
  }
  
  return(best_params)
}

