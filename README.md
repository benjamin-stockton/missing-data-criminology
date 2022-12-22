# Impacts of Incomplete Data on Regression Analysis in Criminology

---

## Overview of the Project

This project investigates the impact of incomplete data on estimation of logistic regression coefficients in the context of race effect estimation for sentencing decisions. We use a simulation study to demonstrate that under different missingness mechanisms it is possible to get results that both over and under estimate the race effect when statistical analysis is performed using complete case analysis, meaning all incomplete cases are dropped from the data set. We further demonstrate that under some missingness mechanisms multiple imputation is able to produce unbiased estimates where complete case analysis cannot.

---

## Contents

This repository is structured as follows:

- Bash scripts to run the simulations are named as "criminology_bash_script_n_###.sh" where '###' can be replaced by the desired sample size. These scripts save simulation outputs to "Sim_Results/simulation_results_Q_#_n_#_m_#.csv" where Q is the simulation iterations, n is the sample size, and m is the number of imputations to perform with multiple imputation. Each of these will take a while to run (a few hours to over a day) so it may be desirable only run some of the simulation settings at a time.
- The simulation script run by the bash scripts is named "criminology_sims.R" and runs a simulation with Q iterations, n sample size, and m imputations where these arguments should be specified in that order from the terminal. See the bash scripts for an example. The script is parallelized and will use $C-1$ cores where $C$ is the total number of cores available. Both the MAR and MNAR mechanisms will be evaluated in a single call.
- The process of determining the missingness mechanism parameters is determined in the file "Data and Missingness Simulation.md".
- The initial data cleaning and preparation with annotations is in "data_cleaning_prep.Rmd". The simulated population is constructed here, while a minimal, reproducible simulation for the population is provided in "simulate_population.R". The "data_cleaning_prep.Rmd" file requires the full PCS data set while the "simulate_population.R" file only requires the $\beta$ coefficients from the PCS logistic regression and data summaries which are saved in the "*.rda" files.
- Some smaller simulated data examples are investigated in the "effect_of_missing_data.Rmd" file to better understand the impacts of different missingness mechanisms.
- "helpers.R" contains all of the essential functions for running the simulations and preparing figures. This source file must be loaded in the majority of the other files.
- "missingness_sims.Rmd" implements a single example of the "criminology_sims.R" simulation with annotations.
- "prepared_figures.R" contains code for recreating figures used in the manuscript.
- "rates_of_missing_information.qmd" implements the multiple imputation study to estimate the rate of missing information in the PCS data and through the rate of missing information determine a range for the number of imputations necessary for a valid analysis.

---

## Requirements

- R (>= 3.4.0)
- R packages:
    - dplyr
    - ggplot2
    - mice
    - magrittr
    - latex2exp
    - RColorBrewer
    - readr
    - cowplot
    - doParallel
    - forcats
    - haven
    - VIM
    - lattice
    - ggthemes
- Quarto

---

## Credits

- Benjamin Stockton; Doctoral Student at the Department of Statistics, University of Connecticut, Storrs, CT
- C. Clare Strange; Assistant Research Professor at the Department of Criminology and Justice Studies, Drexel University, Philadelphia, PA
- Ofer Harel; Associate Dean for Research and Graduate Affairs and Professor of Statistics at the College of Liberal Arts and Sciences, University of Connecticut, Storrs, CT

---
