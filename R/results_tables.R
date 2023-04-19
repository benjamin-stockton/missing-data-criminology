load_sim_results <- function(Q = 225, N = 500, m = 3) {
    for (i in 1:length(m)) {
        fname <- paste0("Sim_Results/simulation_results_Q", Q[i], "_n_", N, "_m_", m[i], ".csv")
        print(fname)
        if (i == 1) {
            sim.res <- read.csv(fname, header = T)
            sim.res$IMPUTATION <- rep(paste0(N[i], " imps"), nrow(sim.res))
        }
        else {
            tmp <- read.csv(fname, header = T)
            tmp$IMPUTATION <- rep(paste0(N[i], " imps"), nrow(tmp))
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

present_tables <- function(Q = 225, N = 500, m = c(3)) {
    
    sim.res <- load_sim_results(Q = Q, N = N, m = m)
    
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
