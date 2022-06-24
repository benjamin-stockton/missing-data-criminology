# sim_functions.R

library(dplyr)
library(doParallel)

generate_data <- function(beta = c(.05, .227, -.2, 0, -.05), M = 1000) {
    CTY <- c(rep(1, M / 12 + 1), rep(2, M / 6), rep(3, M / 2), rep(4, M / 4))
    X <- c(rbinom(M / 12 + 1, 1, .05),
           rbinom(M / 6, 1, .25), 
           rbinom(M / 2, 1, .35),
           rbinom(M / 4, 1, .65))
    C1 <- ifelse(CTY == 1, 1, 0)
    C2 <- ifelse(CTY == 2, 1, 0)
    C4 <- ifelse(CTY == 4, 1, 0)
    XC <- cbind(rep(1, M), X, C1, C2, C4)
    
    EY <- (1 + exp(-XC%*%beta))^(-1)
    Y <- rbinom(M, 1, EY)
    test.df <- data.frame(X = X, Y = Y, C1 = C1, C2 = C2, C4 = C4)
    return(test.df)
}

generate_missing_x <- function(data, miss_type, a = 3, b = 2) {
    M <- nrow(data)
    print(miss_type)
    if (miss_type == "MAR") {
        df.z <- data
        df.z$Z1 <- ifelse(df.z$Y == 0 & df.z$X == 1, 1, 0)
        df.z$Z2 <- ifelse(df.z$Y == 1 & df.z$X == 1, 1, 0)
        
        # CTY predictive of Z1 = (Y == 1 & X == 1)
        fit.z1 <- glm(Z1 ~ C1 + C2 + C4, data = df.z,
                      family = binomial(link = "logit"))
        # CTY predictive of Z2 = (Y == 0 & X == 1)
        fit.z2 <- glm(Z2 ~ C1 + C2 + C4, data = df.z,
                      family = binomial(link = "logit"))
        l1 <- predict(fit.z1, df.z) - 1.7
        p1 <- (1 + exp(-l1))^-1
        l2 <- predict(fit.z2, df.z) - 1.75
        p2 <- (1 + exp(-l2))^-1
        
        p1 <- (1 + exp(-(a + b * ((1 - data$C1) * (1 - data$Y) + (data$C1 * data$Y)))))^-1
        p2 <- (1 + exp(-(a + b * ((1 - data$C1) * data$Y + (data$C1 * (1 - data$Y))))))^-1
        print(table(round(p1, 3)))
        print(table(round(p2, 3)))
    } 
    else {
        p1 <- (1 + exp(-(a + b * ((1 - data$X) * data$Y))))^-1
        p2 <- (1 + exp(-(a + b * ((1 - data$X) * (1 - data$Y)))))^-1
        print(table(round(p1, 3)))
        print(table(round(p2, 3)))
    }
    
    x.miss1 <- rbinom(M, 1, p1)
    x.miss2 <- rbinom(M, 1, p2)
    return(list("COMP" = rep(0, M), "OVER" = x.miss1, "UNDR" = x.miss2))
}

sim_log_reg <- function(data, x.miss, miss_type = "COMP") {
    
    data$X[which(x.miss[[miss_type]] == 1)] <- NA
    p_miss <- mean(is.na(data$X))
    
    fit <- glm(Y ~ X, data = data, family = binomial(link = "logit"))
    # fit <- glm(Y ~ ., data = data, family = binomial(link = "logit"))
    OR <- exp(coef(fit))["X"]
    
    return(c("OR" = OR, "Prop_Miss" = p_miss))
}

full_sim <- function(miss_type = "MAR", beta = c(.05, .227, -.2, 0, -.05),
                     Q = 250, seed = 1234, a = 3, b = 2) {
    set.seed(seed)
    res <- matrix(data = numeric(Q * 12), nrow = Q, ncol = 12)
    levs <- c("comp", "over", "undr")
    colnames(res) <- c(paste0("OR_", levs), 
                       paste0("OR_diff_", levs),
                       paste0("OR_bias_", levs),
                       paste0("PM_", levs))
    
    for (q in 1:Q) {
        if (q %% 50 == 0)
            print(paste0("Iteration: ", q))
        data <- generate_data(beta = beta, M = 1000)
        
        x.miss <- generate_missing_x(data, miss_type, a = a, b = b)
        
        # Complete Data Regression
        res[q, c("OR_comp", "PM_comp")] <- sim_log_reg(data, x.miss, "COMP")
        res[q, c("OR_diff_comp")] <- 0
        res[q, c("OR_bias_comp")] <- res[q, "OR_comp"] - exp(beta[2])
        
        # overerate the race effect under MNAR
        res[q, c("OR_over", "PM_over")] <- sim_log_reg(data, x.miss, "OVER")
        
        res[q, "OR_diff_over"] <- res[q, "OR_over"] - res[q, "OR_comp"]
        res[q, "OR_bias_over"] <- res[q, "OR_over"] - exp(beta[2])
        
        ########
        # Understate the race effect under MNAR
        res[q, c("OR_undr", "PM_undr")] <- sim_log_reg(data, x.miss, "UNDR")
        
        res[q, "OR_diff_undr"] <- res[q, "OR_undr"] - res[q, "OR_comp"]
        res[q, "OR_bias_undr"] <- res[q, "OR_undr"] - exp(beta[2])
    }
    sim.res <- data.frame("ODDS_RATIO" = c(res[,1], res[,2], res[,3]),
                          "OR_DIFF" = c(res[,4], res[,5], res[,6]),
                          "OR_BIAS" = c(res[,7], res[,8], res[,9]),
                          "PROP_MISS" = c(res[,10], res[,11], res[,12]),
                          "MISS_TYPE" = rep(miss_type, 3 * Q),
                          "DIRECTION" = c(rep("COMP", Q),
                                          rep("OVER", Q),
                                          rep("UNDR", Q)))
    
    return(sim.res)
}

result_plots <- function(sim.res) {
    or_plt <- ggplot(sim.res, aes(x = ODDS_RATIO, color = DIRECTION)) +
        geom_density() +
        geom_vline(data = mean_OR, aes(xintercept = odds_ratio_mean, color = DIRECTION)) +
        geom_text(data = mean_OR, 
                  aes(x = odds_ratio_mean, y = 5, label = odds_ratio_mean, color = DIRECTION),
                  nudge_x = .05, angle = 30) +
        labs(title = "Densities of Odds Ratios of the X Effect", 
             x = "Odds Ratio") +
        facet_grid(MISS_TYPE ~ ., scales = "free_y")
    
    pm_plt <- ggplot(sim.res, aes(x = PROP_MISS, color = DIRECTION)) +
        geom_density() +
        geom_vline(data = mean_PM, aes(xintercept = prop_miss_mean, color = DIRECTION)) +
        geom_text(data = mean_PM, 
                  aes(x = prop_miss_mean, y = 45, label = prop_miss_mean, color = DIRECTION),
                  nudge_x = .005, angle = 30) +
        labs(title = "Densities of Prop. of Miss. in X", 
             x = "Missing Proportion") +
        facet_grid(MISS_TYPE ~ ., scales = "free_y")
    
    
    
    diff_plt <- ggplot(sim.res, aes(x = OR_DIFF, color = DIRECTION)) +
        geom_density() +
        geom_vline(data = mean_DIFF, aes(xintercept = odds_ratio_diff_mean, color = DIRECTION)) +
        geom_text(data = mean_DIFF, 
                  aes(x = odds_ratio_diff_mean, y = 2, label = odds_ratio_diff_mean, color = DIRECTION),
                  nudge_x = .01, angle = 30) +
        labs(title = "Densities of Diff in Est. OR of the X Effect", 
             x = TeX("$OR_{miss} - OR_{comp}$")) +
        facet_grid(MISS_TYPE ~ ., scales = "free_y")
    
    bias_plt <- ggplot(sim.res, aes(x = OR_BIAS, color = DIRECTION)) +
        geom_density() +
        geom_vline(data = mean_BIAS, aes(xintercept = bias_mean, color = DIRECTION)) +
        geom_text(data = mean_BIAS, 
                  aes(x = bias_mean, y = 1, label = bias_mean, color = DIRECTION),
                  nudge_x = .05, angle = 30) +
        labs(title = "Densities of Bias in Est. OR of the X Effect", 
             x = TeX("$OR_{miss} - e^\\beta$")) +
        facet_grid(MISS_TYPE ~ ., scales = "free_y")
    
    return(list("OR" = or_plt, "PM" = pm_plt,
                "DIFF" = diff_plt, "BIAS" = bias_plt))
}
