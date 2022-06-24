# sim_functions.R

library(dplyr)
library(doParallel)

generate_data <- function(b = c(.05, .227, -.2, 0, -.05), M = 1000) {
    CTY <- c(rep(1, M / 12 + 1), rep(2, M / 6), rep(3, M / 2), rep(4, M / 4))
    X <- c(rbinom(M / 12 + 1, 1, .05),
           rbinom(M / 6, 1, .25), 
           rbinom(M / 2, 1, .35),
           rbinom(M / 4, 1, .65))
    C1 <- ifelse(CTY == 1, 1, 0)
    C2 <- ifelse(CTY == 2, 1, 0)
    C4 <- ifelse(CTY == 4, 1, 0)
    XC <- cbind(rep(1, M), X, C1, C2, C4)
    
    EY <- (1 + exp(-XC%*%b))^(-1)
    Y <- rbinom(M, 1, EY)
    test.df <- data.frame(X = X, Y = Y, C1 = C1, C2 = C2, C4 = C4)
    return(test.df)
}

generate_missing_x <- function(data, miss_type) {
    M <- nrow(data)
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
        p1 <- predict(fit.z1, df.z, type = "response")
        p2 <- predict(fit.z2, df.z, type = "response")
    } 
    else {
        p1 <- (1 + exp(-(-log(104/4) + 2 * log(5/4) * ((1 - data$X) * data$Y))))^-1
        p2 <- (1 + exp(-(-log(104/4) + 2 * log(5/4) * ((1 - data$X) * (1 - data$Y)))))^-1
    }
    
    x.miss1 <- rbinom(M, 1, p1)
    x.miss2 <- rbinom(M, 1, p2)
    return(list("COMP" = rep(0, M), "EXAG" = x.miss1, "UNDR" = x.miss2))
}

sim_log_reg <- function(data, x.miss, miss_type = "COMP") {
    
    data$X[which(x.miss[[miss_type]] == 1)] <- NA
    p_miss <- mean(is.na(data$X))
    
    fit <- glm(Y ~ ., data = data, family = binomial(link = "logit"))
    OR <- exp(coef(fit))["X"]
    
    return(c("OR" = OR, "Prop_Miss" = p_miss))
}

full_sim <- function(miss_type = "MAR", b = c(.05, .227, -.2, 0, -.05), Q = 250, seed = 1234) {
    set.seed(seed)
    res <- matrix(data = numeric(Q * 12), nrow = Q, ncol = 12)
    levs <- c("comp", "exag", "undr")
    colnames(res) <- c(paste0("OR_", c("comp", "exag", "undr")), 
                       paste0("OR_diff_", c("comp", "exag", "undr")),
                       paste0("OR_bias_", c("comp", "exag", "undr")),
                       paste0("PM_", c("comp", "exag", "undr")))
    
    for (q in 1:Q) {
        if (q %% 50 == 0)
            print(paste0("Iteration: ", q))
        b <- c(.05, 0.227, -.1, 0, .25)
        data <- generate_data(b = b, M = 1000)
        
        x.miss <- generate_missing_x(data, miss_type)
        
        # Complete Data Regression
        res[q, c("OR_comp", "PM_comp")] <- sim_log_reg(data, x.miss, "COMP")
        
        
        # Exagerate the race effect under MNAR
        res[q, c("OR_exag", "PM_exag")] <- sim_log_reg(data, x.miss, "EXAG")
        
        res[q, "OR_diff_exag"] <- res[q, "OR_exag"] - res[q, "OR_comp"]
        res[q, "OR_bias_exag"] <- res[q, "OR_exag"] - exp(b[2])
        
        ########
        # Understate the race effect under MNAR
        res[q, c("OR_undr", "PM_undr")] <- sim_log_reg(data, x.miss, "UNDR")
        
        res[q, "OR_diff_undr"] <- res[q, "OR_undr"] - res[q, "OR_comp"]
        res[q, "OR_bias_undr"] <- res[q, "OR_undr"] - exp(b[2])
    }
    sim.res <- data.frame("ODDS_RATIO" = c(res[,1], res[,2], res[,3]),
                          "OR_DIFF" = c(res[,4], res[,5], res[,6]),
                          "OR_BIAS" = c(res[,7], res[,8], res[,9]),
                          "PROP_MISS" = c(res[,10], res[,11], res[,12]),
                          "MISS_TYPE" = rep(miss_type, 3 * Q),
                          "DIRECTION" = c(rep("COMP", Q),
                                          rep("EXAG", Q),
                                          rep("UNDR", Q)))
    
    return(sim.res)
}