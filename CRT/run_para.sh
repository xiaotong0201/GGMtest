#!/bin/bash

parallel --progress --joblog log.txt -j 100% --results output Rscript Experiment.R {1} {2} ::: 3 7 ::: {1..150}
