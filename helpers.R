# helpers.R

library(dplyr)

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

generate_data <- function(beta = c(.05, .227, -.2, 0, -.05), M = 1000) {
    CTY <- c(rep(1, M / 12 + 1), rep(2, M / 6), rep(3, M / 2), rep(4, M / 4))
    X <- c(rbinom(M / 4, 1, .05),
           rbinom(M / 4, 1, .25), 
           rbinom(M / 4, 1, .35),
           rbinom(M / 4, 1, .65))
    C1 <- ifelse(CTY == 1, 1, 0)
    C2 <- ifelse(CTY == 2, 1, 0)
    C4 <- ifelse(CTY == 4, 1, 0)
    XC <- cbind(rep(1, M), X, C1, C2, C4)
    
    EY <- (1 + exp(-XC%*%beta))^(-1)
    Y <- rbinom(M, 1, EY)
    test.df <- data.frame(X = X, Y = Y, CTY = CTY)
    test.df$CTY <- factor(test.df$CTY)
    return(test.df)
}

missingness_model <- function(data, miss_pars, sim_size = "small") {
    # data: N x p+1 matrix (Y, X)
    # miss_pars: a vector of parameters for the missingness model
    
    data <- dummify_data_matrix(data, sim_size = sim_size)
    
    # print(str(data))
    if (sim_size == "small") {
        V <- data %>% select(-c(X)) %>% as.matrix()
        L <- data %>% select(X) %>% as.matrix()
        W <- L
    }
    else if (sim_size == "full") {
        V <- data %>% select(-c(OFF_RACERBLACK, 
                                OFF_RACERLATINO, 
                                OFF_RACEROTHER)) %>% as.matrix()
        
        L <- data %>% select(c(OFF_RACERBLACK, 
                               OFF_RACERLATINO, 
                               OFF_RACEROTHER)) %>% as.matrix()
        W <- L
    }
    nbeta <- ncol(V); ngamma <- ncol(L); neta <- ncol(W)
    
    alpha <- miss_pars[1]
    beta <- miss_pars[2:(nbeta+1)]
    gamma <- miss_pars[(nbeta+2):(nbeta+ngamma+1)]
    eta <- miss_pars[(nbeta+ngamma+2):(nbeta+ngamma+neta+1)]
    
    # Create probabilities of missingness in Race/X
    logit_miss_prob <- (alpha + V %*% beta + L %*% gamma + W %*% eta)[,1]
    probs <- (1 + exp(-logit_miss_prob))^(-1)
    return(probs)
}

generate_missing_x <- function(data, miss_type, miss_pars, sim_size = "small") {
    M <- nrow(data)
    print(miss_type)
    
    over_pars <- miss_pars[, "OVER"]
    undr_pars <- miss_pars[, "UNDR"]
    
    p1 <- missingness_model(data, miss_pars = over_pars, sim_size = sim_size)
    p2 <- missingness_model(data, miss_pars = undr_pars, sim_size = sim_size)
    print(table(round(p1, 3)))
    print(table(round(p2, 3)))
    
    x.miss1 <- rbinom(M, 1, p1)
    x.miss2 <- rbinom(M, 1, p2)
    return(list("COMP" = rep(0, M), "OVER" = x.miss1, "UNDR" = x.miss2))
}

sim_log_reg <- function(data, x.miss, miss_type = "COMP", sim_size = "small") {
    if (sim_size == "small") {
        data$X[which(x.miss[[miss_type]] == 1)] <- NA
        p_miss <- mean(is.na(data$X))
        
        fit <- glm(Y ~ X, data = data, family = binomial(link = "logit"))
        # fit <- glm(Y ~ ., data = data, family = binomial(link = "logit"))
        OR <- exp(coef(fit))["X"]
        # print(OR)
    }
    else if (sim_size == "full") {
        data$OFF_RACER[which(x.miss[[miss_type]] == 1)] <- NA
        p_miss <- mean(is.na(data$OFF_RACER))
        
        fit <- glm(INCAR ~ ., data = data, family = binomial(link = "logit"))
        # fit <- glm(INCAR ~ OFF_RACER, data = data, family = binomial(link = "logit"))
        OR <- exp(coef(fit))["OFF_RACERBLACK"]
        # print(OR)
    }
    return(c("OR" = OR, "Prop_Miss" = p_miss))
}

full_sim <- function(beta, miss_pars, miss_type = "MAR", 
                     sim_size = "small", pop_data = NULL,
                     Q = 10, replace = FALSE, M = 1000) {
    res <- matrix(data = numeric(Q * 12), nrow = Q, ncol = 12)
    levs <- c("comp", "over", "undr")
    colnames(res) <- c(paste0("OR_", levs),
                       paste0("PM_", levs), 
                       paste0("OR_diff_", levs),
                       paste0("OR_bias_", levs))
    
    for (q in 1:Q) {
        if (q %% 50 == 0)
            print(paste0("Iteration: ", q))
        
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
        x.miss <- generate_missing_x(data, miss_type, sim_size = sim_size,
                                        miss_pars = miss_pars)
        
        # Complete Data Regression
        res[q, c("OR_comp", "PM_comp")] <- sim_log_reg(data, x.miss, "COMP", sim_size = sim_size)
        res[q, c("OR_diff_comp")] <- 0
        res[q, c("OR_bias_comp")] <- res[q, "OR_comp"] - exp(beta.int)
        
        # overerate the race effect under MNAR
        res[q, c("OR_over", "PM_over")] <- sim_log_reg(data, x.miss, "OVER", sim_size = sim_size)
        
        res[q, "OR_diff_over"] <- res[q, "OR_over"] - res[q, "OR_comp"]
        res[q, "OR_bias_over"] <- res[q, "OR_over"] - exp(beta.int)
        
        ########
        # Understate the race effect under MNAR
        res[q, c("OR_undr", "PM_undr")] <- sim_log_reg(data, x.miss, "UNDR", sim_size = sim_size)
        
        res[q, "OR_diff_undr"] <- res[q, "OR_undr"] - res[q, "OR_comp"]
        res[q, "OR_bias_undr"] <- res[q, "OR_undr"] - exp(beta.int)
    }
    sim.res <- data.frame("MISS_TYPE" = rep(miss_type, 3 * Q),
                          "DIRECTION" = c(rep("COMP", Q),
                                          rep("OVER", Q),
                                          rep("UNDR", Q)),
                          "ODDS_RATIO" = c(res[,1], res[,2], res[,3]),
                          "PROP_MISS" = c(res[,4], res[,5], res[,6]),
                          "OR_DIFF" = c(res[,7], res[,8], res[,9]),
                          "OR_BIAS" = c(res[,10], res[,11], res[,12]))
    
    return(sim.res)
}

result_tables <- function(sim.res) {
    mean_OR <- sim.res %>% group_by(MISS_TYPE, DIRECTION) %>%
        summarize(odds_ratio_mean = round(mean(ODDS_RATIO), 3),
                  quant_025 = quantile(ODDS_RATIO, probs = .025),
                  quant_975 = quantile(ODDS_RATIO, probs = .975),
                  or_y_pos = find_label_y_pos(ODDS_RATIO))
    
    mean_PM <- sim.res %>% group_by(MISS_TYPE, DIRECTION) %>%
        summarize(prop_miss_mean = round(mean(PROP_MISS), 3),
                  quant_025 = quantile(PROP_MISS, probs = .025),
                  quant_975 = quantile(PROP_MISS, probs = .975),
                  pm_y_pos = find_label_y_pos(PROP_MISS))
    
    mean_DIFF <- sim.res %>% group_by(MISS_TYPE, DIRECTION) %>%
        summarize(odds_ratio_diff_mean = round(mean(OR_DIFF), 3),
                  quant_025 = quantile(OR_DIFF, probs = .025),
                  quant_975 = quantile(OR_DIFF, probs = .975),
                  diff_y_pos = find_label_y_pos(OR_DIFF))
    
    mean_BIAS <- sim.res %>% group_by(MISS_TYPE, DIRECTION) %>%
        summarize(bias_mean = round(mean(OR_BIAS), 3),
                  quant_025 = quantile(OR_BIAS, probs = .025),
                  quant_975 = quantile(OR_BIAS, probs = .975),
                  bias_y_pos = find_label_y_pos(OR_BIAS))
    return(list("ODDS_RATIO" = mean_OR, "PROP_MISS" = mean_PM, 
                "DIFF" = mean_DIFF, "BIAS" = mean_BIAS))
}

find_label_y_pos <- function(X) {
    d <- density(X)
    return(.5 * max(d$y))
}

result_plots <- function(sim.res) {
    tbls <- result_tables(sim.res)
    mean_OR <- tbls$ODDS_RATIO
    mean_PM <- tbls$PROP_MISS
    mean_DIFF <- tbls$DIFF
    mean_BIAS <- tbls$BIAS
    
    or_plt <- ggplot(sim.res, aes(x = ODDS_RATIO, color = DIRECTION)) +
        geom_density() +
        geom_vline(data = mean_OR, aes(xintercept = odds_ratio_mean, color = DIRECTION)) +
        geom_text(data = mean_OR, 
                  aes(x = odds_ratio_mean, y = or_y_pos,
                      label = odds_ratio_mean, color = DIRECTION),
                  nudge_x = .05, angle = 30) +
        labs(title = "Densities of Odds Ratios of the Race Effect", 
             x = "Odds Ratio") +
        facet_grid(MISS_TYPE ~ ., scales = "free_y")
    
    pm_plt <- ggplot(sim.res, aes(x = PROP_MISS, color = DIRECTION)) +
        geom_density() +
        geom_vline(data = mean_PM, aes(xintercept = prop_miss_mean, color = DIRECTION)) +
        geom_text(data = mean_PM, 
                  aes(x = prop_miss_mean, y = pm_y_pos,
                      label = prop_miss_mean, color = DIRECTION),
                  nudge_x = .005, angle = 30) +
        labs(title = "Densities of Prop. of Miss. in Race", 
             x = "Missing Proportion") +
        facet_grid(MISS_TYPE ~ ., scales = "free_y")
    
    
    diff_plt <- ggplot(sim.res, aes(x = OR_DIFF, color = DIRECTION)) +
        geom_density() +
        geom_vline(data = mean_DIFF, aes(xintercept = odds_ratio_diff_mean, color = DIRECTION)) +
        geom_text(data = mean_DIFF, 
                  aes(x = odds_ratio_diff_mean, y = diff_y_pos,
                      label = odds_ratio_diff_mean, color = DIRECTION),
                  nudge_x = .01, angle = 30) +
        labs(title = "Densities of Diff in Est. OR of the Race Effect", 
             x = TeX("$OR_{miss} - OR_{comp}$")) +
        facet_grid(MISS_TYPE ~ ., scales = "free_y")
    
    bias_plt <- ggplot(sim.res, aes(x = OR_BIAS, color = DIRECTION)) +
        geom_density() +
        geom_vline(data = mean_BIAS, aes(xintercept = bias_mean, color = DIRECTION)) +
        geom_text(data = mean_BIAS, 
                  aes(x = bias_mean, y = bias_y_pos,
                      label = bias_mean, color = DIRECTION),
                  nudge_x = .05, angle = 30) +
        labs(title = "Densities of Bias in Est. OR of the Race Effect", 
             x = TeX("$OR_{miss} - e^\\beta$")) +
        facet_grid(MISS_TYPE ~ ., scales = "free_y")
    
    return(list("OR" = or_plt, "PM" = pm_plt,
                "DIFF" = diff_plt, "BIAS" = bias_plt))
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
