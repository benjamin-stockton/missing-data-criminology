build_miss_par_matrix <- function(beta, gamma, miss_type, sim_size = "small") {
    nvar <- length(beta) + length(gamma)
    miss_pars <- matrix(data = numeric(8 * (nvar+1)), ncol = nvar+1, nrow = 8)
    
    for (i in 1:length(beta)) {
        # column 2 is the response
        miss_pars[,i+1] <- beta[i]
    }
    
    if (sim_size == "small") {
        # gamma
        if (miss_type == "MNAR") {
            miss_pars[,6] <- c(0, rep(gamma[1], 7))
            miss_pars[,7] <- c(0, rep(gamma[2], 7))
            miss_pars[,8] <- c(0, rep(gamma[3], 7))
        }
        else {
            # X
            miss_pars[,6] <- c(0, 0, gamma[1], 0, 3*gamma[1], 0, gamma[1], 0)
            # Z1
            miss_pars[,7] <- c(0, gamma[2], 0, 0, gamma[2], gamma[2], 0, 0)
            # Z1^2
            miss_pars[,8] <- c(0, gamma[3], 0, 0, gamma[3], gamma[3], 0, 0)
            # Z2
            miss_pars[,9] <- c(0, gamma[4], gamma[4], gamma[4], 0, 0, 0, 0)
        }
    }
    else if (sim_size == "full") {
        # gamma
        if (miss_type == "MNAR") {
            for (i in 1:8) {
                miss_pars[, 89+i] <- c(rep(gamma[i], 8))
            }
        }
        else {
            # RACE_BLACK
            miss_pars[,90] <- c(gamma[1], 0, 1.5*gamma[1], 0, 1.5*gamma[1], 0, gamma[1], 0)
            # RACE_LATINO
            miss_pars[,91] <- c(gamma[2], 0, gamma[2], 0, gamma[2], 0, gamma[2], 0)
            # RACE_OTHER
            miss_pars[,92] <- c(gamma[3], 0, gamma[3], 0, gamma[3], 0, gamma[3], 0)
            # AGE
            miss_pars[,93] <- c(0, gamma[4], gamma[4], gamma[4], 0, 0, 0, 0)
            # AGE^2
            miss_pars[,94] <- c(0, gamma[5], gamma[5], gamma[5], 0, 0, 0, 0)
            # RECMIN
            miss_pars[,95] <- c(0, gamma[6], 0, 0, gamma[6], gamma[6], 0, 0)
            # Z1 = I(RACE != WHITE) * (1-INCAR)
            miss_pars[,96] <- c(gamma[7], 0, 1.5*gamma[7], 0, 1.5*gamma[7], 0, gamma[7], 0)
            # Z2 = I(RACE != WHITE) * INCAR
            miss_pars[,97] <- c(gamma[8], 0, 1.5*gamma[8], 0, 1.5*gamma[8], 0, gamma[8], 0)
        }
    }
    
    return(miss_pars)
}