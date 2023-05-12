# criminology_sims_script.R

file.sources = list.files(c("R/sim_funcs"), 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
invisible(capture.output(sapply(file.sources,source,.GlobalEnv)))
invisible(capture.output(source("R/helpers.R")))

# Simulation parameters from command line
args <- c(100, 1000, 3, 0.03)
args <- commandArgs(trailingOnly = T)

if (length(args) != 4) {
    stop("Simulation parameters misspecificed. Should be Q, n, m, p_target")
} else if (length(args) == 4) {
    # Number of iterations to run in each parallel thread
    Q <- as.numeric(args[1]) # Q
    # Sample size
    N <- as.numeric(args[2])
    # Number of imputations for MI
    m <- as.numeric(args[3])
    # Target Proportion of missingness
    p_miss_target <- as.numeric(args[4])
}
print(paste0("Running simulations with ", Q, " iterations, n = ", N, " sample size, and M = ", m, " imputations with p_target = ", p_miss_target))

# Sample with replacement?
replace <- FALSE

# set up cores for parallel computation
cores <- detectCores()[1] - 1
coreQ <- ceiling(Q / cores)

# Load in the simulated data
dat <- read_csv("Data/simulated_data.csv", col_types = "dfdffffffdfdd") %>% as.data.frame()
dat$OFF_RACER <- dat$OFF_RACER %>% forcats::fct_relevel(c("WHITE", "BLACK", "LATINO", "OTHER"))
levels(dat$OFF_RACER)

# load in the coefficients from the CCA logistic regression on the full data set
# named 'beta'
load("Data/full_data_cca_log_reg_coef.rda")

beta["OFF_RACERBLACK"]
# sort(beta, decreasing = T)[1:5]
# sort(beta)[1:10]


##------------------------------------------------------------------------------
# Test simulation
print("Test Simulation: 1 Iteration of MAR")

# Missingness model parameters for MNAR - OVER
miss_pars_over <- build_miss_par_matrix(
    beta = log(c(10, rep(1, 87))),
    gamma = log(c(25, rep(1, 5), 1,1)),
    miss_type = "MAR", sim_size = "full")
# print(miss_pars_over[, c(1:4, 90:97)])

# Missingness model parameters for MNAR - UNDR
miss_pars_undr <- build_miss_par_matrix(
    beta = log(c(.1, rep(1, 87))),
    gamma = log(c(1, rep(1, 5), 25, .1)),
    miss_type = "MAR", sim_size = "full")

fs <- full_sim(miss_type = "MAR", sim_size = "full", pop_data = dat,
               p_miss_target = p_miss_target,
               beta = beta, N = N, m = m, miss_pars_over = miss_pars_over,
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
        invisible(capture.output(sapply(file.sources,source,.GlobalEnv)))
        invisible(capture.output(source("R/helpers.R")))
        full_sim(miss_type = "MNAR", sim_size = "full", pop_data = dat,
                 p_miss_target = p_miss_target,
                 beta = beta, N = N, m = m, miss_pars_over = miss_pars_mnar_over, 
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
        invisible(capture.output(sapply(file.sources,source,.GlobalEnv)))
        invisible(capture.output(source("R/helpers.R")))
        full_sim(miss_type = "MAR", sim_size = "full", pop_data = dat,
                 p_miss_target = p_miss_target,
                 beta = beta, N = N, m = m, miss_pars_over = miss_pars_mar_over, 
                 miss_pars_undr = miss_pars_mar_undr, Q = coreQ, replace = replace)
    }
})[3]
ptime
stopCluster(cl)
mar

sim.res <- rbind(mnar, mar)

write_csv(sim.res, paste0("Sim_Results/simulation_results_Q", Q, "_n_", N, 
                          "_m_", m, "_p_miss_", p_miss_target, ".csv"))
