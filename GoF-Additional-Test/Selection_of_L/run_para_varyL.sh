#!/bin/bash

parallel --progress --joblog log.txt -j 100% --results output Rscript Experiment_varyL.R {1} {2} ::: {5..5} ::: {1..400}
