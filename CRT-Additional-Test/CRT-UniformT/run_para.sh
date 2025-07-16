#!/bin/bash

parallel --progress --joblog log.txt -j 100% --results output Rscript Experiment.R {1} {2} ::: 1 5 ::: {1..400}
