
full_sim <- function(beta, miss_pars_over, miss_pars_undr,
                     p_miss_target = 0.01, miss_type = "MAR", 
                     sim_size = "small", pop_data = NULL,
                     Q = 1, replace = FALSE, N = 100, m = 3) {
    res <- matrix(data = numeric(Q * 12), nrow = Q, ncol = 12)
    levs <- c("COMP", "OVER", "UNDR")
    lbls <- character(0)
    for (l in levs) {
        lbls <- c(
              lbls, 
              paste0(
                  c("OR_", "PM_", "OR_diff_", "OR_bias_"), 
                  l)
              )
    }
    colnames(res) <- lbls
    
    mi_res <- matrix(data = numeric(Q * 8), nrow = Q, ncol = 8)
    
    colnames(mi_res) <- lbls[5:12]
    
    for (q in 1:Q) {
        if (q %% 50 == 0)
            print(paste0("Iteration: ", q))
        # Simulation Step
        tmp <- simulate_analysis(beta = beta,
                                 miss_pars_over = miss_pars_over,
                                 miss_pars_undr = miss_pars_undr, 
                                 sim_size = sim_size,
                                 pop_data = pop_data,
                                 replace = replace, 
                                 N = N,
                                 m = m,
                                 levs = levs)
        res[q,] <- tmp$cc
        mi_res[q,] <- tmp$mi
    }
    sim.res <- data.frame("MISS_TYPE" = rep(miss_type, 3 * Q),
                          "DIRECTION" = c(rep("COMP", Q),
                                          rep("OVER", Q),
                                          rep("UNDR", Q)),
                          "ODDS_RATIO" = c(res[,1], res[,5], res[,9]),
                          "PROP_MISS" = c(res[,2], res[,6], res[,10]),
                          "OR_DIFF" = c(res[,3], res[,7], res[,11]),
                          "OR_BIAS" = c(res[,4], res[,8], res[,12]),
                          "ANALYSIS" = rep("CCA", 3*Q))
    
    mi.res <- data.frame("MISS_TYPE" = rep(miss_type, 2 * Q),
                         "DIRECTION" = c(rep("OVER", Q),
                                         rep("UNDR", Q)),
                         "ODDS_RATIO" = c(mi_res[,1], mi_res[,5]),
                         "PROP_MISS" = c(mi_res[,2], mi_res[,6]),
                         "OR_DIFF" = c(mi_res[,3], mi_res[,7]),
                         "OR_BIAS" = c(mi_res[,4], mi_res[,8]),
                         "ANALYSIS" = rep("MI", 2*Q))
    sim.res <- rbind(sim.res, mi.res)
    
    return(sim.res)
}

simulate_analysis <- function(beta, miss_pars_over, miss_pars_undr,
                              p_miss_target = 0.01, 
                              sim_size = "small", pop_data = NULL,
                              replace = FALSE, N = 100, m = 3,
                              levs = c("COMP", "OVER", "UNDR")) {
    
    # Generate the Data or Sample from the Population Data
    if (sim_size == "full") {
        if (is.null(pop_data)) {
            print("Oops, no population data provided!")
            break
        }
        data <- sample_from_pop(pop_dat = pop_data,
                                N = N,
                                replace = replace)
        
        beta.int <- beta["OFF_RACERBLACK"]
    }
    else if (sim_size == "small") {
        data <- generate_data(beta = beta, N = N)
        beta.int <- beta["X"]
    }
    
    # Create the missing data
    dat.list <- generate_missing_x(data, p_miss_target = p_miss_target,
                                   sim_size = sim_size,
                                   miss_pars_over = miss_pars_over,
                                   miss_pars_undr = miss_pars_undr)
    
    # Run Logistic Regression and get results for each setting:
    #   COMP: Complete data
    #   OVER: Overstate the effect
    #   UNDR: Understate the effect
    res <- numeric(0)
    mi_res <- numeric(0)
    
    for (l in levs) {
        lr_res <- sim_log_reg(dat.list[[l]], sim_size = sim_size)
        if (l == "COMP") {
            cc_or <- lr_res[1]
        }
        cc_pm <- lr_res[2]
        res <- c(res, 
                 lr_res[1], # OR
                 lr_res[2], # PM
                 lr_res[1] - cc_or,         # OR_Diff
                 lr_res[1] - exp(beta.int)) # OR_Bias
        if (l != "COMP") {
            lr_res <- mi_log_reg(dat.list[[l]], sim_size = sim_size, m = m)
            mi_res <- c(mi_res, 
                        lr_res, # OR
                        cc_pm,     # PM
                        lr_res - cc_or,         # OR_Diff
                        lr_res - exp(beta.int)) # OR_Bias
        }
    }
    return(list("cc" = res, "mi" = mi_res))
}