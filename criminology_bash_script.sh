#!/usr/bin/env bash

echo "Running big sims"
Rscript criminology_sims.R 105 10000 25
Rscript criminology_sims.R 105 25000 25
Rscript criminology_sims.R 105 50000 25

Rscript criminology_sims.R 60 100000 3

