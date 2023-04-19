#!/usr/bin/env bash

echo "Running n = 500 sims"
Rscript R/criminology_sims.R 225 500 3
Rscript R/criminology_sims.R 225 500 5
Rscript R/criminology_sims.R 225 500 8
Rscript R/criminology_sims.R 225 500 25
Rscript R/criminology_sims.R 225 500 50
Rscript R/criminology_sims.R 225 500 100