library(ggplot2)
library(cowplot)
library(dplyr)

set.seed(987)
M <- 1000
CTY <- c(rep(1, M / 4), rep(2, M / 4), rep(3, M / 4), rep(4, M / 4))
X <- c(rbinom(M / 4, 1, .05),
       rbinom(M / 4, 1, .15), 
       rbinom(M / 4, 1, .25),
       rbinom(M / 4, 1, .5))
b <- 0.227

EY <- exp(.05+b*X) / (1 + exp(.05+b*X))
table(EY)
Y <- numeric(M)
for (i in 1:M) {
    Y[i] <- rbinom(1, 1, EY[i])
}
test.df <- data.frame(X = X, Y = Y, CTY = factor(CTY))

fit3 <- glm(Y ~ X + CTY, data = test.df, family = binomial(link = "logit"))
# summary(fit3)
print("Coef")
coef(fit3)
print("OR")
exp(coef(fit3))[2]

generate_missing_x <- function(data, miss_type, dir, prop_miss1, prop_miss2) {
    M <- nrow(data)
    x.miss <- rep(0, M)
    if (dir == "EXAG" & miss_type == "MNAR") {
        for (i in 1:M) {
            if (data$X[i] == 1 & data$Y[i] == 0) {
                x.miss[i] <- rbinom(1, 1, prop_miss1)
            }
        }
    }
    else if (dir == "UNDR" & miss_type == "MNAR") {
        for (i in 1:M) {
            if (data$X[i] == 1 & data$Y[i] == 1) {
                x.miss[i] <- rbinom(1, 1, prop_miss1)
            }
        }
    }
    else if (dir == "EXAG" & miss_type == "MAR") {
        for (i in 1:M) {
            if ((data$Y[i] == 0 & (data$CTY[i] == 3 | data$CTY[i] == 4)) |
                (data$Y[i] == 1 & data$CTY[i] == 1)) {
                x.miss[i] <- rbinom(1, 1, prop_miss1)
            }
        }
    }
    else {
        for (i in 1:M) {
            if ((data$Y[i] == 0 & data$CTY[i] == 1) |
                (data$Y[i] == 1 & data$CTY[i] == 4)) {
                x.miss[i] <- rbinom(1, 1, prop_miss1)
            }
        }
    }
    return(x.miss)
}

sim_log_reg_w_miss <- function(prop_miss1 = .28, prop_miss2 = .2,
                               miss_type = "MAR", Q = 250, seed = 1234) {
    set.seed(seed)
    OR_true <- numeric(Q)
    OR_exag <- numeric(Q)
    OR_undr <- numeric(Q)
    OR_diff_true <- rep(0, Q)
    OR_diff_exag <- numeric(Q)
    OR_diff_undr <- numeric(Q)
    OR_bias_true <- numeric(Q)
    OR_bias_exag <- numeric(Q)
    OR_bias_undr <- numeric(Q)
    p_miss_true <- rep(0, Q)
    p_miss_exag <- numeric(Q)
    p_miss_undr <- numeric(Q)
    
    for (q in 1:Q) {
        M <- 1000
        CTY <- c(rep(1, M / 4), rep(2, M / 4), rep(3, M / 4), rep(4, M / 4))
        X <- c(rbinom(M / 4, 1, .05),
               rbinom(M / 4, 1, .15), 
               rbinom(M / 4, 1, .25),
               rbinom(M / 4, 1, .5))
        b <- 0.227
        
        EY <- exp(.05+b*X) / (1 + exp(.05+b*X))
        table(EY)
        Y <- numeric(M)
        for (i in 1:M) {
            Y[i] <- rbinom(1, 1, EY[i])
        }
        data <- data.frame(X = X, Y = Y, CTY = factor(CTY))
        fit <- glm(Y ~ X, data = data, family = binomial(link = "logit"))
        OR_true[q] <- exp(coef(fit))[2]
        OR_bias_true[q] <- OR_true[q] - exp(b)
        rm(fit)
        
        # Exagerate the race effect under MNAR
        if (q %% 50 == 0)
            print(paste0("Iteration: ", q))
        df1 <- data
        df2 <- data
        
        x.miss <- generate_missing_x(df1, miss_type, dir = "EXAG",
                                     prop_miss1 = prop_miss1,
                                     prop_miss2 = prop_miss2)
        df1$X[which(x.miss == 1)] <- NA
        p_miss_exag[q] <- mean(is.na(df1$X))
        
        fit <- glm(Y ~ X, data = df1, family = binomial(link = "logit"))
        OR_exag[q] <- exp(coef(fit))[2]
        OR_diff_exag[q] <- OR_exag[q] - OR_true[q]
        OR_bias_exag[q] <- OR_exag[q] - exp(b)
        rm(df1)
        rm(fit)
        
        ########
        # Understate the race effect under MNAR
        
        x.miss <- generate_missing_x(df2, miss_type, dir = "UNDR",
                                     prop_miss1 = prop_miss1,
                                     prop_miss2 = prop_miss2)
        df2$X[which(x.miss == 1)] <- NA
        p_miss_undr[q] <- mean(is.na(df2$X))
        
        fit <- glm(Y ~ X, data = df2, family = binomial(link = "logit"))
        OR_undr[q] <- exp(coef(fit))[2]
        OR_diff_undr[q] <- OR_undr[q] - OR_true[q]
        OR_bias_undr[q] <- OR_undr[q] - exp(b)
        rm(fit)
        rm(df2)
    }
    sim.res <- data.frame("ODDS_RATIO" = c(OR_true, OR_exag, OR_undr),
                          "OR_DIFF" = c(OR_diff_true, OR_diff_exag, OR_diff_undr),
                          "OR_BIAS" = c(OR_bias_true, OR_bias_exag, OR_bias_undr),
                          "PROP_MISS" = c(p_miss_true, p_miss_exag, p_miss_undr),
                          "MISS_TYPE" = rep(miss_type, 3 * Q),
                          "DIRECTION" = c(rep("TRUE", Q),
                                          rep("EXAG", Q),
                                          rep("UNDR", Q)))
    
    return(sim.res)
}


### MNAR

# Want to get roughly 3% missing in X in both cases.

# Exagerate the X effect
# Missingness in X | X = 1 & Y = 0, P(X is miss | X = 1, Y = 0) = .28

# Understate the X effect
# Missingness in X | X = 1 & Y = 1, P(X is miss | X = 1, Y = 1) = .2


mnar <- sim_log_reg_w_miss(prop_miss1 = .3, prop_miss2 = .225, miss_type = "MNAR")

##########################################

# MAR

# Want to get roughly 3% missing in X in both cases.

# Exagerate the X effect
# Missingness in X | CTY = 3,4 & Y = 0, P(X is miss | Y = 1) = .28

# Understate the X effect
# Missingness in X | CTY = 3,4 & Y = 1, P(X is miss | Y = 0) = .2

mar <- sim_log_reg_w_miss(prop_miss1 = .12, prop_miss2 = .2, miss_type = "MAR")

#########################
# Results

sim.res <- rbind(mnar, mar)

mean_OR <- sim.res %>% group_by(DIRECTION, MISS_TYPE) %>%
    summarize(mean_val = round(mean(ODDS_RATIO), 3),
              quant_025 = quantile(ODDS_RATIO, probs = .025),
              quant_975 = quantile(ODDS_RATIO, probs = .975))
mean_OR

mean_PM <- sim.res %>% group_by(DIRECTION, MISS_TYPE) %>%
    summarize(mean_val = round(mean(PROP_MISS), 3),
              quant_025 = quantile(PROP_MISS, probs = .025),
              quant_975 = quantile(PROP_MISS, probs = .975))
mean_PM

mean_DIFF <- sim.res %>% group_by(DIRECTION, MISS_TYPE) %>%
    summarize(mean_val = round(mean(OR_DIFF), 3),
              quant_025 = quantile(OR_DIFF, probs = .025),
              quant_975 = quantile(OR_DIFF, probs = .975))
mean_DIFF

mean_BIAS <- sim.res %>% group_by(DIRECTION, MISS_TYPE) %>%
    summarize(mean_val = round(mean(OR_BIAS), 3),
              quant_025 = quantile(OR_BIAS, probs = .025),
              quant_975 = quantile(OR_BIAS, probs = .975))
mean_BIAS

or_plt <- ggplot(sim.res, aes(x = ODDS_RATIO, color = DIRECTION)) +
    geom_density() +
    geom_vline(data = mean_OR, aes(xintercept = mean_val, color = DIRECTION)) +
    geom_text(data = mean_OR, 
              aes(x = mean_val, y = 6, label = mean_val, color = DIRECTION),
              nudge_x = .075, angle = 30) +
    labs(title = "Densities of Odds Ratios of the X Effect", 
         x = "Odds Ratio") +
    facet_grid(MISS_TYPE ~ DIRECTION)

pm_plt <- ggplot(sim.res, aes(x = PROP_MISS, color = DIRECTION)) +
    geom_density() +
    geom_vline(data = mean_PM, aes(xintercept = mean_val, color = DIRECTION)) +
    geom_text(data = mean_PM, 
              aes(x = mean_val, y = 45, label = mean_val, color = DIRECTION),
              nudge_x = .005, angle = 30) +
    labs(title = "Densities of Prop. of Miss. in X", 
         x = "Missing Proportion") +
    facet_grid(MISS_TYPE ~ DIRECTION)

diff_plt <- ggplot(sim.res, aes(x = OR_DIFF, color = DIRECTION)) +
    geom_density() +
    geom_vline(data = mean_DIFF, aes(xintercept = mean_val, color = DIRECTION)) +
    geom_text(data = mean_DIFF, 
              aes(x = mean_val, y = 5, label = mean_val, color = DIRECTION),
              nudge_x = .005, angle = 30) +
    labs(title = "Densities of Diff in Est. OR of the X Effect", 
         x = "OR_miss - OR_comp") +
    facet_grid(MISS_TYPE ~ DIRECTION)

bias_plt <- ggplot(sim.res, aes(x = OR_BIAS, color = DIRECTION)) +
    geom_density() +
    geom_vline(data = mean_BIAS, aes(xintercept = mean_val, color = DIRECTION)) +
    geom_text(data = mean_BIAS, 
              aes(x = mean_val, y = 1, label = mean_val, color = DIRECTION),
              nudge_x = .05, angle = 30) +
    labs(title = "Densities of Bias in Est. OR of the X Effect", 
         x = "OR_miss - e^\beta") +
    facet_grid(MISS_TYPE ~ DIRECTION)

cowplot::plot_grid(or_plt, pm_plt, nrow = 2)
cowplot::plot_grid(diff_plt, bias_plt, nrow = 2)

or_bias_diff <- ggplot(sim.res, aes(x = OR_BIAS, y = OR_DIFF, color = DIRECTION)) +
    geom_point() +
    facet_grid(MISS_TYPE ~ DIRECTION)

cowplot::plot_grid(bias_plt, or_bias_diff, nrow = 2)
