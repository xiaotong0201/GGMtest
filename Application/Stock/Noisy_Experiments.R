

## add noise experiments

load("DataForExp.RData")

sub_Y=Y 
setT=setT2

library(parallel)

noise_level= c(0.05, 0.10, 0.15, 0.20, 0.25)


num_repeats <- 400


run_noisy = function(epo) {
  source('FunctionsStock.R')
  set.seed(epo)
 
  
  exp_Y <- sub_Y + rnorm(length(sub_Y), sd = noise_sd)
  
  
  return(para_run_one(mat_X, exp_Y, setT, G1))
}

# run_noisy(1)   # test before running

############# Parallel begin #####
noisy_results = list()

# Detect the number of cores
no_cores <-
  detectCores()   # leave one core free for system processes

# Create a cluster
cl <- makeCluster(no_cores)

for (noise_sd in noise_level  ) {
  cat("Processing noise sd:", noise_sd, "\n")
  
  clusterExport(cl,
                c("run_noisy", "num_repeats",
                  'G1', 'mat_X', 'setT', 'sub_Y', 'noise_sd'))
  
  new_mc_results <- parLapply(cl, 1:num_repeats , run_noisy)
  noisy_results[[as.character(noise_sd)]] <- new_mc_results
}
stopCluster(cl)



save(noisy_results, num_repeats,noise_level, file = 'Stock-noise-DIA-significant.RData')


sum(unlist(lapply(noisy_results, function(x)sapply(x, `[[`, 'use.time')))) / 4


noisy_pvalues_table=simplify2array(lapply(noisy_results, function(x)
  sapply(x, `[[`, 'pvalues')))

sig.lv=0.05
apply(noisy_pvalues_table,c(3,1), function(v){mean(v<sig.lv) } )
