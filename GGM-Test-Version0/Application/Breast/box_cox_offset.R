library(nortest)

# Transformation function
transform_offset <- function(x, lambda, offset=10) {
  if(lambda==0){
    log(x+offset)
  }else{
    (((x+offset))^lambda -1)/lambda
  }
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
      transformed_x <- transform_offset(x, lambda, offset=gamma)
      p_value <- test_normality(transformed_x)
      
      if (p_value > best_p_value) {
        best_p_value <- p_value
        best_params <- c(lambda, gamma)
      }
    }
  }
  
  return(best_params)
}


