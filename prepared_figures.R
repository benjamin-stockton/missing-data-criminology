source("helpers.R")
library(magrittr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(forcats)

#####################################################
# Loading the results
#####################################################

sim.res.500 <- load_sim_results(Q = 225, M = 500, m = c(3, 5, 8, 25))
head(sim.res.500)

sim.res.500 %<>% mutate(IMPUTATION = ifelse(ANALYSIS == "MI", levels(IMPUTATION)[IMPUTATION], ""),
                        ANALYSIS2 = ifelse(ANALYSIS == "MI", paste0(ANALYSIS, " - ", IMPUTATION), levels(ANALYSIS)[ANALYSIS]),
                       SampleSize = "N = 500")

head(sim.res.500)
sim.res.500 %>% filter(ANALYSIS == "MI") %>% head()
sim.res.500 %>% filter(ANALYSIS == "MI") %>% tail()

sim.res.1000 <- load_sim_results(Q = 225, M = 1000, m = c(3, 5, 8, 25))
sim.res.1000 %<>% mutate(IMPUTATION = ifelse(ANALYSIS == "MI", levels(IMPUTATION)[IMPUTATION], ""),
                        ANALYSIS2 = ifelse(ANALYSIS == "MI", paste0(ANALYSIS, " - ", IMPUTATION), levels(ANALYSIS)[ANALYSIS]),
                        SampleSize = "N = 1000")

sim.res.2500 <- load_sim_results(Q = 225, M = 2500, m = c(3, 5, 8))
sim.res.2500 %<>% mutate(IMPUTATION = ifelse(ANALYSIS == "MI", levels(IMPUTATION)[IMPUTATION], ""),
                        ANALYSIS2 = ifelse(ANALYSIS == "MI", paste0(ANALYSIS, " - ", IMPUTATION), levels(ANALYSIS)[ANALYSIS]),
                        SampleSize = "N = 2500")

sim.res.25002 <- load_sim_results(Q = 105, M = 2500, m = 25)
sim.res.25002 %<>% mutate(IMPUTATION = ifelse(ANALYSIS == "MI", levels(IMPUTATION)[IMPUTATION], ""),
                         ANALYSIS2 = ifelse(ANALYSIS == "MI", paste0(ANALYSIS, " - ", IMPUTATION), levels(ANALYSIS)[ANALYSIS]),
                         SampleSize = "N = 2500")

sim.res.25000 <- load_sim_results(Q = 105, M = 25000, m = 3)
sim.res.25000 %<>% mutate(IMPUTATION = ifelse(ANALYSIS == "MI", levels(IMPUTATION)[IMPUTATION], ""),
                        ANALYSIS2 = ifelse(ANALYSIS == "MI", paste0(ANALYSIS, " - ", IMPUTATION), levels(ANALYSIS)[ANALYSIS]),
                        SampleSize = "N = 25000")

#####################################################
# Formatting Results for Plotting
#####################################################

sim.res <- bind_rows(sim.res.500, sim.res.1000, sim.res.2500, sim.res.25002)

table(sim.res$SampleSize)

sim.res$DIRECTION %<>% forcats::fct_relevel(c("UNDR", "No Missing Data", "OVER")) %>% forcats::fct_recode(Under = "UNDR", Over = "OVER")
sim.res$SampleSize %<>% forcats::fct_relevel(paste0("N = ", c(500, 1000, 2500, 25000)))
sim.res$ANALYSIS2 %<>% forcats::fct_relevel(c("COMP", "CCA", paste0("MI - ", c(3, 5, 8, 25), " imps")))

med_comp <- median(sim.res[sim.res$ANALYSIS == "COMP", "ODDS_RATIO"])

#####################################################
# Plotting Results
#####################################################

# Odds Ratio
ggplot(sim.res, aes(DIRECTION, ODDS_RATIO, fill = SampleSize)) +
    geom_hline(yintercept = med_comp, color = "gray") +
    geom_boxplot(alpha = .75) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_x", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Odds Ratio")

ggsave2("Plots/sim_res_all_sample_sizes_odds_ratio.png", dpi = 600, width = 8, height = 6, units = "in")

# ggplot(sim.res, aes(DIRECTION, ODDS_RATIO, fill = ANALYSIS2)) +
#     geom_hline(yintercept = med_comp) +
#     geom_boxplot() +
#     facet_grid(MISS_TYPE ~ SampleSize) +
#     theme_bw() +
#     scale_fill_grey()

# Plot bias
ggplot(sim.res, aes(DIRECTION, OR_BIAS, fill = SampleSize)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_x", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Stat. Bias of the Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Statistical Bias")

ggsave2("Plots/sim_res_all_sample_sizes_or_bias.png", dpi = 600, width = 8, height = 6, units = "in")

# Plot missing cases proportions

ggplot(sim.res.nc, aes(DIRECTION, PROP_MISS, fill = SampleSize)) +
    geom_boxplot(alpha = 0.75) +
    facet_wrap(MISS_TYPE~.) +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Proportion of Incomplete Cases", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Proportion of Incomplete Cases")

ggsave2("Plots/sim_res_all_sample_sizes_prop_miss.png", dpi = 600, width = 8, height = 6, units = "in")

# Filter out Complete Data Analyses
sim.res.nc <- sim.res %>% filter(ANALYSIS != "COMP")

ggplot(sim.res.nc, aes(DIRECTION, OR_DIFF, fill = SampleSize)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_x", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Difference in Race Effect Estimates from CCA/MI to the Complete Data Estimate", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Difference in the OR")

ggsave2("Plots/sim_res_all_sample_sizes_or_diff.png", dpi = 600, width = 8, height = 6, units = "in")



# N = 25000 Plots and table

sim.res.25000$DIRECTION %<>% forcats::fct_relevel(c("UNDR", "No Missing Data", "OVER")) %>% forcats::fct_recode(Under = "UNDR", Over = "OVER")
sim.res.25000 %>% group_by(DIRECTION, ANALYSIS2, MISS_TYPE) %>% summarize(
    mean_OR = mean(ODDS_RATIO),
    sd_OR = sd(ODDS_RATIO),
    mean_Diff = mean(OR_DIFF),
    sd_Diff = sd(OR_DIFF),
    mean_bias = mean(OR_BIAS),
    sd_bias = sd(OR_BIAS)
)

ggplot(sim.res.25000, aes(DIRECTION, ODDS_RATIO)) +
    geom_hline(yintercept = med_comp, color = "gray") +
    geom_boxplot(alpha = .75) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_x", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Odds Ratio")

ggplot(sim.res.25000, aes(DIRECTION, OR_DIFF)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_x", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Difference in Race Effect Estimates from CCA/MI to the Complete Data Estimate", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Difference in the OR")

ggplot(sim.res.25000, aes(DIRECTION, OR_BIAS)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_x", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Stat. Bias of the Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Statistical Bias")

