# helpers.R

library(dplyr, warn.conflicts = FALSE)

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
library(mice)

sample_from_pop <- function(pop_dat, M = 1000, replace = FALSE) {
    s <- sample(1:nrow(pop_dat), size = M, replace = replace)
    samp.dat <- pop_dat[s,]
    return(samp.dat)
}

dummify_data_matrix <- function(data, sim_size = "small") {
    if (sim_size == "small") {
        data %<>% mutate(as.data.frame(model.matrix(~ 0 + CTY, data = cur_data()))) %>% 
            select(-c(CTY, CTY1))
        return(data)
        
    }
    else if (sim_size == "full") {
        data %<>%
            mutate(as.data.frame(model.matrix(~ 0 + OFF_RACER + PRS + CRIMETYPE + YEAR + COUNTY, data = cur_data()))) %>%
            select(-c(COUNTY, OFF_RACER, PRS, YEAR, CRIMETYPE, OFF_RACERWHITE)) 
        data$TRIAL <- ifelse(data$TRIAL == "Yes", 1, 0)
        data$MALE <- ifelse(data$MALE == "Male", 1, 0)
        data$RECMIN <- ifelse(data$RECMIN == "Yes", 1, 0)
        
        return(data)
    }
}

generate_data <- function(beta = c(.05, .227, 0, -.01, .9, -.2, 0, -.05), M = 1000) {
    CTY <- sample(1:4, size = M, replace = TRUE, prob = c(.05, .25, .25, .45))
    X <- c(rbinom(M / 4, 1, .05),
           rbinom(M / 4, 1, .25), 
           rbinom(M / 4, 1, .35),
           rbinom(M / 4, 1, .65))
    C1 <- ifelse(CTY == 1, 1, 0)
    C2 <- ifelse(CTY == 2, 1, 0)
    C4 <- ifelse(CTY == 4, 1, 0)
    Z1 <- rnorm(M)
    Z1Q <- Z1^2
    Z2 <- rbinom(M, 1, .45)
    XC <- cbind(rep(1, M), X, Z1, Z1Q, Z2, C1, C2, C4)
    
    EY <- (1 + exp(-XC%*%beta))^(-1)
    Y <- rbinom(M, 1, EY)
    test.df <- data.frame(X = X, Y = Y, Z1 = Z1, Z1Q = Z1Q, Z2 = Z2, CTY = CTY)
    test.df$CTY <- factor(test.df$CTY)
    return(test.df)
}

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

calc_pattern_probs <- function(V, L, beta, gamma, EV, EL) {
    pat.probs <- c(834548, 26206, 2078, 34, 1465, 92, 1, 1) / 864422
    alpha <- numeric(8)
    prob.mat <- matrix(numeric(8 * nrow(V)), ncol = 8)
    
    for (i in 1:8) {
        alpha[i] <- -log(1 / pat.probs[i] - 1) - (EV %*% beta[i,] + EL %*% gamma[i,])[,1]
        lp_i <- alpha[i] + (V %*% beta[i,] + L %*% gamma[i,])[,1]
        p_i <- exp(lp_i)
        prob.mat[,i] <- p_i
    }
    prob.mat[,1] <- rep(1, nrow(prob.mat))
    r.sum <- apply(prob.mat, 1, sum)
    
    nmat <- t(apply(matrix(1:nrow(prob.mat), ncol = 1), 1, function(r) {
        return(prob.mat[r,] / r.sum[r])
    }))
    # nmat[,1] <- sapply(p1, function(p) {return(ifelse(p > 0, p, 0))})
    
    return(nmat)
}

missingness_model <- function(data, miss_pars, sim_size = "small") {
    # data: N x p+1 matrix (Y, X)
    # miss_pars: a matrix of parameters for the missingness model
    
    data <- dummify_data_matrix(data, sim_size = sim_size)
    
    # print(str(data))
    if (sim_size == "small") {
        V <- data %>% select(-c(X, Z1, Z1Q, Z2)) %>% as.matrix()
        L <- data %>% select(X, Z1, Z1Q, Z2) %>% as.matrix()
    }
    else if (sim_size == "full") {
        V <- data %>% select(-c(OFF_RACERBLACK, 
                                OFF_RACERLATINO, 
                                OFF_RACEROTHER, 
                                DOSAGE, DOSAGEQ,
                                RECMIN)) %>% as.matrix()
        
        L <- data %>% select(c(INCAR, OFF_RACERBLACK, 
                               OFF_RACERLATINO, 
                               OFF_RACEROTHER, 
                               DOSAGE, DOSAGEQ,
                               RECMIN)) %>% 
            mutate(Z1 = OFF_RACERBLACK * (1 - INCAR),
                   Z2 = OFF_RACERBLACK * INCAR) %>% select(-INCAR) %>% as.matrix()
    }
    # print(colnames(V[,1:6]))
    # print(colnames(L))
    nbeta <- ncol(V); ngamma <- ncol(L)
    
    alpha <- miss_pars[,1]
    beta <- miss_pars[,2:(nbeta+1)]
    gamma <- miss_pars[,(nbeta+2):(nbeta+ngamma+1)]
    
    EV <- apply(V, 2, mean); EL <- apply(L, 2, mean)
    prob.mat <- calc_pattern_probs(V, L, beta, gamma, EV, EL)
    
    # Create probabilities of missingness in Race/X
    assigned_pattern <- apply(matrix(1:nrow(data), ncol = 1), 1,
                              function(r) {
                                  return(sample(1:8, 1, prob = prob.mat[r,]))
                              })
    # print(prop.table(table(assigned_pattern)))
    return(assigned_pattern)
}

assign_missing <- function(data, miss_pattern, sim_size = "small") {
    df <- data
    
    if (sim_size == "small") {
        df[which(miss_pattern == 2), "X"] <- NA
        df[which(miss_pattern == 3), "Z1"] <- NA
        df[which(miss_pattern == 4), c("X", "Z1")] <- NA
        df[which(miss_pattern == 5), "Z2"] <- NA
        df[which(miss_pattern == 6), c("Z2", "X")] <- NA
        df[which(miss_pattern == 7), c("Z2", "Z1")] <- NA
        df[which(miss_pattern == 8), c("Z1", "Z2", "X")] <- NA
        
        df$Z1Q <- df$Z1^2
    } 
    else if (sim_size == "full") {
        df[which(miss_pattern == 2), "OFF_RACER"] <- NA
        df[which(miss_pattern == 3), "RECMIN"] <- NA
        df[which(miss_pattern == 4), c("OFF_RACER", "RECMIN")] <- NA
        df[which(miss_pattern == 5), "DOSAGE"] <- NA
        df[which(miss_pattern == 6), c("DOSAGE", "OFF_RACER")] <- NA
        df[which(miss_pattern == 7), c("DOSAGE", "RECMIN")] <- NA
        df[which(miss_pattern == 8), c("DOSAGE", "RECMIN", "OFF_RACER")] <- NA
        
        # Induce missingness in dosage^2
        df$DOSAGEQ <- df$DOSAGE^2
    }
    
    return(df)
}

generate_missing_x <- function(data, miss_pars_over, miss_pars_undr, sim_size = "small") {
    p1 <- missingness_model(data, miss_pars = miss_pars_over, sim_size = sim_size)
    p2 <- missingness_model(data, miss_pars = miss_pars_undr, sim_size = sim_size)
    data.over <- assign_missing(data, p1, sim_size = sim_size)
    data.under <- assign_missing(data, p2, sim_size = sim_size)
    
    return(list("COMP" = data, "OVER" = data.over, "UNDR" = data.under))
}

sim_log_reg <- function(data, sim_size = "small") {
    df <- data %>% filter(complete.cases(.))
    # if (sim_size == "full") {
    #     Z1 <- ifelse(df$OFF_RACER == "BLACK" & df$INCAR == 1, 1, 0)
    #     Z2 <- ifelse(df$OFF_RACER == "WHITE" & df$INCAR == 1, 1, 0)
    #     Z3 <- ifelse(df$OFF_RACER == "BLACK" & df$INCAR == 0, 1, 0)
    #     Z4 <- ifelse(df$OFF_RACER == "WHITE" & df$INCAR == 0, 1, 0)
    #     
    #     # print(paste0("P(RACE == BLACK & INCAR == 1) = ", mean(Z1)))
    #     # p_miss <- (mean(Z1) / mean(Z3)) / (mean(Z2) / mean(Z4))
    # }
    # else {
    #     p_miss <- 1 - nrow(df) / nrow(data)
    # }
    p_miss <- 1 - nrow(df) / nrow(data)
    
    
    if (sim_size == "small") {
        fit <- glm(Y ~ ., data = df, family = binomial(link = "logit"))
        OR <- exp(coef(fit))["X"]
    }
    else if (sim_size == "full") {
        # print(prop.table(table(df$INCAR, df$OFF_RACER), margin = 2))
        fit <- glm(INCAR ~ ., data = df, family = binomial(link = "logit"))
        OR <- exp(coef(fit))["OFF_RACERBLACK"]
    }
    return(c("OR" = OR, "Prop_Miss" = p_miss))
}

mi_log_reg <- function(data, sim_size = "small", m = 3) {
    # Create the multpily imputed data sets
    imp <- mice(data, m = m)
    # Stack the data sets together
    # imp_tot <- complete(imp, "long", inc = TRUE)
    # Run the analyses
    if (sim_size == "small") {
        fitm <- with(imp, glm(Y ~ X + Z1 + Z1Q + Z2 + CTY, family = binomial(link = "logit")))
        pool_fit <- pool(fitm)
        x_index <- which(pool_fit$pooled[,1] == "X")
        OR <- exp(pool_fit$pooled[x_index, "estimate"])
    }
    else if (sim_size == "full") {
        fitm <- with(imp, glm(INCAR ~ CRIMETYPE + OGS + OGSQ + RECMIN + TRIAL + PRS + 
                              MALE + DOSAGE + DOSAGEQ + OFF_RACER + COUNTY + YEAR, family = binomial(link = "logit")))
        pool_fit <- pool(fitm)
        race_index <- which(pool_fit$pooled[,1] == "OFF_RACERBLACK")
        OR <- exp(pool_fit$pooled[race_index, "estimate"])
    }
    return(OR)
}

full_sim <- function(beta, miss_pars_over, miss_pars_undr, miss_type = "MAR", 
                     sim_size = "small", pop_data = NULL,
                     Q = 1, replace = FALSE, M = 100, m = 3) {
    res <- matrix(data = numeric(Q * 12), nrow = Q, ncol = 12)
    levs <- c("COMP", "OVER", "UNDR")
    colnames(res) <- c(paste0("OR_", levs),
                       paste0("PM_", levs), 
                       paste0("OR_diff_", levs),
                       paste0("OR_bias_", levs))
    mi_res <- matrix(data = numeric(Q * 8), nrow = Q, ncol = 8)
    levs2 <- c("OVER", "UNDR")
    colnames(mi_res) <- c(paste0("OR_", levs2),
                          paste0("PM_", levs2),
                       paste0("OR_diff_", levs2),
                       paste0("OR_bias_", levs2))
    
    for (q in 1:Q) {
        if (q %% 50 == 0)
            print(paste0("Iteration: ", q))
        
        # Generate the Data or Sample from the Population Data
        if (sim_size == "full") {
            if (is.null(pop_data)) {
                print("Oops, no population data provided!")
                break
            }
            data <- sample_from_pop(pop_dat = pop_data, M = M,
                                    replace = replace)
            
            beta.int <- beta["OFF_RACERBLACK"]
        }
        else if (sim_size == "small") {
            data <- generate_data(beta = beta, M = M)
            beta.int <- beta["X"]
        }
        
        # Create the missing data
        dat.list <- generate_missing_x(data, sim_size = sim_size,
                                        miss_pars_over = miss_pars_over,
                                       miss_pars_undr = miss_pars_undr)
        
        # Run Logistic Regression and get results for each setting:
        #   COMP: Complete data
        #   OVER: Overstate the effect
        #   UNDR: Understate the effect
        for (l in levs) {
            lr_res <- sim_log_reg(dat.list[[l]], sim_size = sim_size)
            res[q, c(paste0("OR_", l), paste0("PM_", l))] <- lr_res
            res[q, paste0("OR_diff_", l)] <- res[q, paste0("OR_", l)] - res[q, "OR_COMP"]
            res[q, paste0("OR_bias_", l)] <- res[q, paste0("OR_", l)] - exp(beta.int)
            if (l != "COMP") {
                lr_res <- mi_log_reg(dat.list[[l]], sim_size = sim_size, m = m)
                mi_res[q, paste0("OR_", l)] <- lr_res
                mi_res[q, paste0("PM_", l)] <- res[q, paste0("PM_", l)]
                mi_res[q, paste0("OR_diff_", l)] <- mi_res[q, paste0("OR_", l)] - res[q, "OR_COMP"]
                mi_res[q, paste0("OR_bias_", l)] <- mi_res[q, paste0("OR_", l)] - exp(beta.int)
            }
        }
    }
    sim.res <- data.frame("MISS_TYPE" = rep(miss_type, 3 * Q),
                          "DIRECTION" = c(rep("COMP", Q),
                                          rep("OVER", Q),
                                          rep("UNDR", Q)),
                          "ODDS_RATIO" = c(res[,1], res[,2], res[,3]),
                          "PROP_MISS" = c(res[,4], res[,5], res[,6]),
                          "OR_DIFF" = c(res[,7], res[,8], res[,9]),
                          "OR_BIAS" = c(res[,10], res[,11], res[,12]),
                          "ANALYSIS" = rep("CCA", 3*Q))
    mi.res <- data.frame("MISS_TYPE" = rep(miss_type, 2 * Q),
                         "DIRECTION" = c(rep("OVER", Q),
                                         rep("UNDR", Q)),
                         "ODDS_RATIO" = c(mi_res[,1], mi_res[,2]),
                         "PROP_MISS" = c(mi_res[,3], mi_res[,4]),
                         "OR_DIFF" = c(mi_res[,5], mi_res[,6]),
                         "OR_BIAS" = c(mi_res[,7], mi_res[,8]),
                         "ANALYSIS" = rep("MI", 2*Q))
    sim.res <- rbind(sim.res, mi.res)
    
    return(sim.res)
}

load_sim_results <- function(Q = 225, M = 500, m = 3) {
    for (i in 1:length(m)) {
        fname <- paste0("Sim_Results/simulation_results_Q", Q, "_n_", M, "_m_", m[i], ".csv")
        if (i == 1) {
            sim.res <- read.csv(fname, header = T)
            sim.res$IMPUTATION <- rep(paste0(m[i], " imps"), nrow(sim.res))
        }
        else {
            tmp <- read.csv(fname, header = T)
            tmp$IMPUTATION <- rep(paste0(m[i], " imps"), nrow(tmp))
            sim.res <- rbind(sim.res, tmp)
        }
    }
    sim.res[which(sim.res$DIRECTION == "COMP"), "ANALYSIS"] <- "COMP"
    sim.res[which(sim.res$DIRECTION == "COMP"), "DIRECTION"] <- "No Missing Data"
    sim.res$DIRECTION %<>% fct_relevel(c("OVER", "UNDR", "No Missing Data"))
    sim.res$ANALYSIS %<>% fct_relevel(c("CCA", "COMP", "MI"))
    sim.res$IMPUTATION %<>% fct_relevel(paste0(m, " imps"))
    
    return(sim.res)
}

result_tables <- function(sim.res) {
    options(dplyr.summarise.inform = FALSE)
    mean_OR <- sim.res %>% group_by(IMPUTATION, MISS_TYPE, DIRECTION, ANALYSIS) %>%
        summarise(odds_ratio_mean = round(mean(ODDS_RATIO), 3),
                  quant_025 = quantile(ODDS_RATIO, probs = .025),
                  quant_975 = quantile(ODDS_RATIO, probs = .975),
                  or_y_pos = find_label_y_pos(ODDS_RATIO))
    
    mean_PM <- sim.res %>% filter(ANALYSIS != "COMP") %>%
        group_by(MISS_TYPE, DIRECTION) %>%
        summarise(prop_miss_mean = round(mean(PROP_MISS), 3),
                  quant_025 = round(quantile(PROP_MISS, probs = .025), 3),
                  quant_975 = round(quantile(PROP_MISS, probs = .975), 3),
                  pm_y_pos = find_label_y_pos(PROP_MISS))
    
    mean_DIFF <- sim.res %>% group_by(IMPUTATION, MISS_TYPE, DIRECTION, ANALYSIS) %>%
        summarise(odds_ratio_diff_mean = round(mean(OR_DIFF), 3),
                  quant_025 = quantile(OR_DIFF, probs = .025),
                  quant_975 = quantile(OR_DIFF, probs = .975),
                  diff_y_pos = find_label_y_pos(OR_DIFF))
    
    mean_BIAS <- sim.res %>% group_by(IMPUTATION, MISS_TYPE, DIRECTION, ANALYSIS) %>%
        summarise(bias_mean = round(mean(OR_BIAS), 3),
                  quant_025 = quantile(OR_BIAS, probs = .025),
                  quant_975 = quantile(OR_BIAS, probs = .975),
                  bias_y_pos = find_label_y_pos(OR_BIAS))
    return(list("ODDS_RATIO" = mean_OR, "PROP_MISS" = mean_PM, 
                "DIFF" = mean_DIFF, "BIAS" = mean_BIAS))
}

present_tables <- function(Q = 225, M = 500, m = c(3)) {
    
    sim.res <- load_sim_results(Q = Q, M = M, m = m)
    
    tbls <- result_tables(sim.res)
    
    cnames <- c("Num of Imp", "Miss Type", "Direction", "Analysis", "Odds Ratio", "0.025 Quantile", "0.975 Quantile")
    stat_names <- c("Odds Ratio", "Prop. of Inc. Cases", "OR Difference", "Bias")
    for (i in 1:4) {
        if (i == 2) {
            tmp <- tbls[[i]][,1:5]
            cnames[5] <- stat_names[i]
            colnames(tmp) <- cnames[c(2:3, 5:7)]
        }
        else {
            tmp <- tbls[[i]][,1:7]
            cnames[5] <- stat_names[i]
            colnames(tmp) <- cnames
        }
        
        tbls[[i]] <- tmp
    }
    return(tbls)
}

find_label_y_pos <- function(X) {
    d <- density(X)
    return(.5 * max(d$y))
}

result_plots <- function(sim.res, Q = 225, M = 1000, horizontal = T) {
    # options(dplyr.summarise.inform = FALSE)
    # tbls <- result_tables(sim.res)
    # mean_OR <- tbls$ODDS_RATIO
    # mean_PM <- tbls$PROP_MISS
    # mean_DIFF <- tbls$DIFF
    # mean_BIAS <- tbls$BIAS
    
    or_plt <- my_gg_box(sim.res, X = "ODDS_RATIO", 
                        main = "Dist. of Odds Ratios of the Race Effect",
                        xlab = "Odds Ratio", Q = Q, M = M, horizontal = horizontal)
    
    bias_plt <- my_gg_box(sim.res, X = "OR_BIAS",
                          main = "Densities of Bias in Est. OR of the Race Effect", 
                          xlab = TeX("$OR_{miss} - e^\\beta$"), Q = Q, M = M, horizontal = horizontal)
    
    # Removing the complete data estimates from the table since they aren't needed for the other plots
    sim.res <- sim.res[which(sim.res$ANALYSIS != "COMP"),]
    
    pm_plt <- my_gg_box(sim.res, X = "PROP_MISS", 
                         main = "Dist. of Prop. of Missing Cases",
                        xlab = "Missing Proportion", Q = Q, M = M, horizontal = horizontal)
    
    diff_plt <- my_gg_box(sim.res, X = "OR_DIFF", 
                          main = "Densities of Diff in Est. OR of the Race Effect",
                          xlab = TeX("$OR_{miss} - OR_{comp}$"), Q = Q, M = M, horizontal = horizontal)
    
    return(list("OR" = or_plt, "PM" = pm_plt,
                "DIFF" = diff_plt, "BIAS" = bias_plt))
}

present_figs <- function(M = 500, Q = 225, m = c(3,5,8), horizontal = horizontal, ...) {
    
    sim.res <- load_sim_results(Q = Q, M = M, m = m)
    
    plt.list <- result_plots(sim.res, Q = Q, M = M, horizontal = horizontal)
    
    return(plt.list)
}

my_gg_bar <- function(dat, X, ...) {
    p1 <- ggplot(dat, aes_string(x = X), ) +
        geom_bar(aes(y = prop.table(stat(count))), 
                 stat = "count") + 
        geom_text(aes(y = prop.table(stat(count)), 
                      label = scales::percent(prop.table(stat(count)))),
                  stat = 'count',
                  position = position_dodge(.9),
                  vjust = -0.5,
                  size = 3) +
        scale_y_continuous(labels = scales::percent) + 
        labs(y = "Percent")
    return(p1)
}

my_gg_box <- function(dat, X, main = "", xlab = "", Q = 225, M = 500, horizontal = T, ...) {
    p1 <- ggplot(dat, aes_string(x = X)) +
        geom_vline(xintercept = 0, color = "gray") +
        geom_boxplot(aes(y = IMPUTATION, fill = DIRECTION)) +
        labs(title = main,
             subtitle = paste0(Q, " iterations, n = ", M),
             x = xlab) +
        facet_grid(MISS_TYPE ~ ANALYSIS, scales = "free_y") +
        scale_fill_grey()
    
    if (!horizontal) {
        p1 <- p1 + coord_flip()
    }
    return(p1)
}

# Density function is broken
# my_gg_dens <- function(dat, X, main = "", x = "", ...) {
#     p1 <- ggplot(dat, aes_string(x = X)) +
#         geom_density(aes(y = ANALYSIS, color = DIRECTION)) +
#         geom_vline(data = mean_X, aes(xintercept = odds_ratio_diff_mean, color = DIRECTION)) +
#         geom_text(data = mean_X,
#                   aes(x = odds_ratio_diff_mean, y = diff_y_pos,
#                       label = odds_ratio_diff_mean, color = DIRECTION),
#                   nudge_x = .01, angle = 30) +
#         labs(title = "Densities of Diff in Est. OR of the Race Effect",
#              subtitle = paste0(Q, " iterations, n = ", M, ", M = ", m),
#              x = TeX("$OR_{miss} - OR_{comp}$")) +
#         facet_grid(MISS_TYPE ~ ., scales = "free_y") +
#         scale_color_colorblind()
#     return(p1)
# }
