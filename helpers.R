# helpers.R

library(dplyr)

generate_missing_x <- function(data, miss_type, dir, prop_miss) {
    race_by_county <- data %>% group_by(COUNTY) %>%
        summarise(prop_white = mean(OFF_RACER == "WHITE", na.rm = T),
                  prop_black = mean(OFF_RACER == "BLACK", na.rm = T),
                  prop_latino = mean(OFF_RACER == "LATINO", na.rm = T),
                  prop_other = mean(OFF_RACER == "OTHER", na.rm = T)) %>%
        as.data.frame()
    
    cty_set_1 <- race_by_county %>% filter(prop_white < .55) %>%
        arrange(desc(prop_black)) %>% select(COUNTY)
    cty_set_2 <- race_by_county %>% filter(prop_white > .75) %>%
        arrange(desc(prop_white)) %>% select(COUNTY)
    cty_set_1 <- cty_set_1[,1]
    cty_set_2 <- cty_set_2[,1]
    
    M <- nrow(data)
    x.miss <- rep(0, M)
    if (dir == "OVER" & miss_type == "MNAR") {
        for (i in 1:M) {
            if (data$OFF_RACER[i] != "WHITE" & data$INCAR[i] == 0) {
                x.miss[i] <- rbinom(1, 1, prop_miss)
            }
        }
    }
    else if (dir == "UNDR" & miss_type == "MNAR") {
        for (i in 1:M) {
            if (data$OFF_RACER[i] != "WHITE" & data$INCAR[i] == 1) {
                x.miss[i] <- rbinom(1, 1, prop_miss)
            }
        }
    }
    else if (dir == "OVER" & miss_type == "MAR") {
        for (i in 1:M) {
            if ((data$INCAR[i] == 0 & (data$COUNTY[i] %in% cty_set_1)) |
                (data$INCAR[i] == 1 & data$COUNTY[i] %in% cty_set_2)) {
                x.miss[i] <- rbinom(1, 1, prop_miss)
            }
        }
    }
    else {
        for (i in 1:M) {
            if ((data$INCAR[i] == 0 & data$COUNTY[i] %in% cty_set_2) |
                (data$INCAR[i] == 1 & data$COUNTY[i] %in% cty_set_1)) {
                x.miss[i] <- rbinom(1, 1, prop_miss)
            }
        }
    }
    return(x.miss)
}

sim_log_reg_w_miss <- function(data, prop_miss1 = .28, prop_miss2 = .2,
                               miss_type = "MAR", Q = 100, seed = 1234) {
    set.seed(seed)
    OR_true <- numeric(Q)
    OR_OVER <- numeric(Q)
    OR_UNDR <- numeric(Q)
    OR_diff_true <- rep(0, Q)
    OR_diff_OVER <- numeric(Q)
    OR_diff_UNDR <- numeric(Q)
    OR_bias_true <- numeric(Q)
    OR_bias_OVER <- numeric(Q)
    OR_bias_UNDR <- numeric(Q)
    p_miss_true <- rep(0, Q)
    p_miss_OVER <- numeric(Q)
    p_miss_UNDR <- numeric(Q)
    
    # Known beta
    b <- 0.227
    
    for (q in 1:Q) {
        
        fit <- glm(INCAR ~ ., data = data, family = binomial(link = "logit"))
        # print(summary(fit))
        OR_true[q] <- exp(coef(fit))["OFF_RACERBLACK"]
        # print(log(OR_true[q]))
        OR_bias_true[q] <- OR_true[q] - exp(b)
        rm(fit)
        
        # OVERerate the race effect under MNAR
        if (q %% 25 == 0)
            print(paste0("Iteration: ", q))
        df1 <- data
        df2 <- data
        
        x.miss <- generate_missing_x(df1, miss_type, dir = "OVER",
                                     prop_miss = prop_miss1)
        df1$OFF_RACER[which(x.miss == 1)] <- NA
        p_miss_OVER[q] <- mean(is.na(df1$OFF_RACER))
        
        fit <- glm(INCAR ~ ., data = df1, family = binomial(link = "logit"))
        # print(summary(fit))
        OR_OVER[q] <- exp(coef(fit))["OFF_RACERBLACK"]
        # print(log(OR_OVER[q]))
        OR_diff_OVER[q] <- OR_OVER[q] - OR_true[q]
        OR_bias_OVER[q] <- OR_OVER[q] - exp(b)
        rm(df1)
        rm(fit)
        
        ########
        # Understate the race effect under MNAR
        
        x.miss <- generate_missing_x(df2, miss_type, dir = "UNDR",
                                     prop_miss = prop_miss2)
        df2$OFF_RACER[which(x.miss == 1)] <- NA
        p_miss_UNDR[q] <- mean(is.na(df2$OFF_RACER))
        
        fit <- glm(INCAR ~ ., data = df2, family = binomial(link = "logit"))
        # print(summary(fit))
        OR_UNDR[q] <- exp(coef(fit))["OFF_RACERBLACK"]
        # print(log(OR_UNDR[q]))
        OR_diff_UNDR[q] <- OR_UNDR[q] - OR_true[q]
        OR_bias_UNDR[q] <- OR_UNDR[q] - exp(b)
        rm(fit)
        rm(df2)
    }
    sim.res <- data.frame("ODDS_RATIO" = c(OR_true, OR_OVER, OR_UNDR),
                          "OR_DIFF" = c(OR_diff_true, OR_diff_OVER, OR_diff_UNDR),
                          "OR_BIAS" = c(OR_bias_true, OR_bias_OVER, OR_bias_UNDR),
                          "PROP_MISS" = c(p_miss_true, p_miss_OVER, p_miss_UNDR),
                          "MISS_TYPE" = rep(miss_type, 3 * Q),
                          "DIRECTION" = c(rep("TRUE", Q),
                                          rep("OVER", Q),
                                          rep("UNDR", Q)))
    
    return(sim.res)
}

