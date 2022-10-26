# criminology_sims_script.R

# Load packages
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
library(magrittr)
library(ggplot2)
library(latex2exp)
library(RColorBrewer)
library(readr)
library(cowplot)
library(doParallel)
library(forcats)
source("helpers.R")

# Simulation parameters from command line
args <- commandArgs(trailingOnly = T)

if (length(args) != 3) {
    stop("Simulation parameters misspecificed. Should be Q, n, m")
} else if (length(args) == 3) {
    # Number of iterations to run in each parallel thread
    Q <- as.numeric(args[1]) # Q
    # Sample size
    M <- as.numeric(args[2])
    # Number of imputations for MI
    m <- as.numeric(args[3])
}
print(paste0("Running simulations with ", Q, " iterations, n = ", M, " sample size, and M = ", m, " imputations"))

# Sample with replacement?
replace <- FALSE

# set up cores for parallel computation
cores <- detectCores()[1] - 1
coreQ <- ceiling(Q / cores)

# Load in the simulated data
dat <- read_csv("Data/simulated_data.csv", col_types = "dfdffffffdfdd") %>% as.data.frame()
dat$OFF_RACER %<>% fct_relevel(c("WHITE", "BLACK", "LATINO", "OTHER"))
levels(dat$OFF_RACER)

# load in the coefficients from the CCA logistic regression on the full data set
load("full_data_cca_log_reg_coef.rda")

beta["OFF_RACERBLACK"]
# sort(beta, decreasing = T)[1:5]
# sort(beta)[1:10]


##------------------------------------------------------------------------------
# Test simulation
print("Test Simulation: 1 Iteration of MNAR")

# Missingness model parameters for MNAR - OVER
miss_pars_over <- build_miss_par_matrix(
    beta = log(c(1, rep(1, 87))),
    gamma = log(c(1, rep(1, 5), .001, 25)),
    miss_type = "MNAR", sim_size = "full")
# print(miss_pars_over[, c(1:4, 90:97)])

# Missingness model parameters for MNAR - UNDR
miss_pars_undr <- build_miss_par_matrix(
    beta = log(c(1, rep(1, 87))),
    gamma = log(c(1, rep(1, 5), 25, .1)),
    miss_type = "MNAR", sim_size = "full")

fs <- full_sim(miss_type = "MNAR", sim_size = "full", pop_data = dat,
               beta = beta, M = M, m = m, miss_pars_over = miss_pars_over,
               miss_pars_undr = miss_pars_undr, Q = 1, replace = replace)
fs

##------------------------------------------------------------------------------
# MNAR simulation
print(paste0("MNAR Simulation: ", Q, " iterations"))

# Missingness model parameters for MNAR - OVER
miss_pars_mnar_over <- build_miss_par_matrix(
    beta = log(c(10, rep(1, 87))),
    gamma = log(c(.1, rep(1, 5), .001, 100)),
    miss_type = "MNAR", sim_size = "full")

# Missingness model parameters for MNAR - UNDR
miss_pars_mnar_undr <- build_miss_par_matrix(
    beta = log(c(5, rep(1, 87))),
    gamma = log(c(5, rep(1, 5), 10, 1)),
    miss_type = "MNAR", sim_size = "full")

cl <- makeCluster(cores)
registerDoParallel(cl)
ptime <- system.time({
    mnar <- foreach(i = 1:cores, .combine = rbind) %dopar% {
        source("helpers.R")
        full_sim(miss_type = "MNAR", sim_size = "full", pop_data = dat,
                 beta = beta, M = M, m = m, miss_pars_over = miss_pars_mnar_over, 
                 miss_pars_undr = miss_pars_mnar_undr, Q = coreQ, replace = replace)
    }
})[3]
ptime
stopCluster(cl)

##------------------------------------------------------------------------------
# MAR simulation
print(paste0("MAR Simulation: ", Q, " iterations"))

# Missingness model parameters for MAR - OVER
miss_pars_mar_over <- build_miss_par_matrix(
    beta = log(c(10, rep(1, 87))),
    gamma = log(c(.1, rep(1, 5), .00001, 1)),
    miss_type = "MAR", sim_size = "full")

# Missingness model parameters for MAR - UNDR
miss_pars_mar_undr <- build_miss_par_matrix(
    beta = log(c(1, rep(1, 87))),
    gamma = log(c(1, rep(1, 5), 30, .1)),
    miss_type = "MAR", sim_size = "full")

cl <- makeCluster(cores)
registerDoParallel(cl)
ptime <- system.time({
    mar <- foreach(i = 1:cores, .combine = rbind) %dopar% {
        source("helpers.R")
        full_sim(miss_type = "MAR", sim_size = "full", pop_data = dat,
                 beta = beta, M = M, m = m, miss_pars_over = miss_pars_mar_over, 
                 miss_pars_undr = miss_pars_mar_undr, Q = coreQ, replace = replace)
    }
})[3]
ptime
stopCluster(cl)

sim.res <- rbind(mnar, mar)

write_csv(sim.res, paste0("Sim_Results/simulation_results_Q", Q, "_n_", M, "_m_", m, ".csv"))
