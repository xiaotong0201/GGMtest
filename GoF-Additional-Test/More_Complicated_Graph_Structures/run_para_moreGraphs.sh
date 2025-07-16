#!/bin/bash

# Define parameters
MODELS=4   # Number of models
REPS=400   # Number of replications
CORES=15   # Number of parallel jobs

# Generate and run commands dynamically using GNU parallel
parallel --progress --joblog log.txt --results output -j$CORES --env Rscript --colsep ' ' Rscript run_experiments.R {1} {2} ::: $(seq 4 $MODELS) ::: $(seq 1 $REPS)
