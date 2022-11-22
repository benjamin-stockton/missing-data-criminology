source("helpers.R")
library(magrittr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(forcats)

#####################################################
# Loading the results
#####################################################

sim_res_list <- list()

sim_settings <- read.csv("Sim_Results/sim_settings.csv", header = T)

sim_set2 <- sim_settings# %>% filter(n < 30000)
Ns <- unique(sim_set2$n)


for (n in Ns) {
    M <- sim_set2[which(sim_set2$n == n), "M"]; Q <- sim_set2[which(sim_set2$n == n), "simQ"]
    
    tmp_str <- paste0("res_", n)
    print(tmp_str)
    print(Q)
    tmp <- load_sim_results(Q = Q, N = n, M = M)

    tmp %<>% mutate(IMPUTATION = ifelse(ANALYSIS == "MI", levels(IMPUTATION)[IMPUTATION], ""),
                    ANALYSIS2 = ifelse(ANALYSIS == "MI", paste0(ANALYSIS, " - ", IMPUTATION), levels(ANALYSIS)[ANALYSIS]),
                    SampleSize = paste0("N = ", n))

    sim_res_list[[tmp_str]] <- tmp
}

#####################################################
# Formatting Results for Plotting
#####################################################

# sim.res <- bind_rows(sim.res.500, sim.res.1000, sim.res.2500, sim.res.25002)
sim.res <- purrr::reduce(sim_res_list, bind_rows)
table(sim.res$SampleSize)

sim.res$DIRECTION %<>% forcats::fct_relevel(c("UNDR", "No Missing Data", "OVER")) %>% forcats::fct_recode(Under = "UNDR", Over = "OVER")
sim.res$SampleSize %<>% forcats::fct_relevel(paste0("N = ", c(500, 1000, 2500, 5000, 10000, 25000, 50000)))
sim.res$ANALYSIS2 %<>% forcats::fct_relevel(c("COMP", "CCA", paste0("MI - ", c(3, 5, 8, 25), " imps")))

med_comp <- median(sim.res[sim.res$ANALYSIS == "COMP", "ODDS_RATIO"])

smy.res <- sim.res %>% group_by(SampleSize, MISS_TYPE, DIRECTION, ANALYSIS2) %>% summarize(
    mean_OR = mean(ODDS_RATIO),
    sd_OR = sd(ODDS_RATIO),
    mean_Diff = mean(OR_DIFF),
    sd_Diff = sd(OR_DIFF),
    min_Diff = min(OR_DIFF),
    max_Diff = max(OR_DIFF),
    mean_bias = mean(OR_BIAS),
    sd_bias = sd(OR_BIAS)
)
smy.res %>% print(n = 66)

sim.res %>% group_by(MISS_TYPE, SampleSize, DIRECTION, ANALYSIS2) %>% summarize(
    mean_OR = mean(ODDS_RATIO),
    sd_OR = sd(ODDS_RATIO),
    mean_Diff = mean(OR_DIFF),
    sd_Diff = sd(OR_DIFF),
    mean_bias = mean(OR_BIAS),
    sd_bias = sd(OR_BIAS)
)%>% mutate(
    "Avg OR (sd)" = paste0(round(mean_OR, 3), " (", round(sd_OR, 3), ")"),
    "Avg Diff (sd)" = paste0(round(mean_Diff, 3), " (", round(sd_Diff, 3), ")"),
    "Avg Bias (sd)" = paste0(round(mean_bias, 3), " (", round(sd_bias, 3), ")")
) %>% 
    select(MISS_TYPE, SampleSize, DIRECTION, ANALYSIS2, 'Avg OR (sd)', 'Avg Diff (sd)', 'Avg Bias (sd)') %>% write.csv(file = "Sim_Results/summarized_sim_results.csv", row.names = F)

sim.res %>% filter (DIRECTION %in% c("Under", "Over")) %>%
    group_by(SampleSize, MISS_TYPE, DIRECTION) %>% summarize(
    "Avg. Prop of Incomplete Cases" = mean(PROP_MISS),
    "0.025 Quantile" = quantile(PROP_MISS, c(0.025)),
    "0.975 Quantile" = quantile(PROP_MISS, c(0.975))
) %>% write.csv(file = "Sim_Results/summarized_prop_miss_results.csv", row.names = F)

sim.res %>% filter(DIRECTION == "Under", ANALYSIS2 == "CCA") %>%
    group_by(SampleSize, MISS_TYPE, DIRECTION, ANALYSIS2) %>% summarize(
    mean_OR = mean(ODDS_RATIO),
    sd_OR = sd(ODDS_RATIO),
    "Avg OR (sd)" = paste0(round(mean_OR, 3), " (", round(sd_OR, 3), ")"),
    "0.025 Quantile" = round(quantile(ODDS_RATIO, c(0.025)), 3),
    "0.975 Quantile" = round(quantile(ODDS_RATIO, c(0.975)), 3),
    "95 CI" = paste0("(", `0.025 Quantile`, ",", `0.975 Quantile`, ")")
    
) %>% select(SampleSize, MISS_TYPE, DIRECTION, ANALYSIS2, 'Avg OR (sd)', '95 CI') %>% write.csv(file = "Sim_Results/summarized_cca_under_results.csv", row.names = F)

sim.res %>% filter(DIRECTION == "Under", ANALYSIS2 == "CCA", MISS_TYPE == "MAR") %>%
    group_by(SampleSize, MISS_TYPE, DIRECTION, ANALYSIS2) %>% summarize(
        mean_OR = mean(ODDS_RATIO),
        sd_OR = sd(ODDS_RATIO),
        mean_Diff = mean(OR_DIFF),
        sd_Diff = sd(OR_DIFF),
        "Avg OR (sd)" = paste0(round(mean_OR, 3), " (", round(sd_OR, 3), ")"),
        "0.025 Quantile" = round(quantile(ODDS_RATIO, c(0.025)), 3),
        "0.975 Quantile" = round(quantile(ODDS_RATIO, c(0.975)), 3),
        "95 CI" = paste0("(", `0.025 Quantile`, ",", `0.975 Quantile`, ")")
        
    )

sim.res %>% filter(DIRECTION == "Under", ANALYSIS2 == "CCA") %>%
    group_by(MISS_TYPE, DIRECTION, ANALYSIS2) %>% summarize(
        mean_OR = mean(ODDS_RATIO),
        sd_OR = sd(ODDS_RATIO),
        mean_Diff = mean(OR_DIFF),
        sd_Diff = sd(OR_DIFF),
        "Avg OR (sd)" = paste0(round(mean_OR, 3), " (", round(sd_OR, 3), ")"),
        "0.025 Quantile" = round(quantile(ODDS_RATIO, c(0.025)), 3),
        "0.975 Quantile" = round(quantile(ODDS_RATIO, c(0.975)), 3),
        "95 CI" = paste0("(", `0.025 Quantile`, ",", `0.975 Quantile`, ")")
        
    )

#####################################################
# Plotting Results
#####################################################

# Small Sample Size Plots
plot_filepath <- "Plots/sim_res_small_sample_sizes_"
small_smpl.res <- sim.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps") & SampleSize %in% paste0("N = ", c(500, 1000, 2500, 5000)))

# Odds Ratio
ggplot(small_smpl.res, aes(DIRECTION, ODDS_RATIO, fill = SampleSize)) +
    geom_hline(yintercept = med_comp, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1, outlier.shape = 1) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_x", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Odds Ratio")

ggsave2(paste0(plot_filepath, "odds_ratio.png"), dpi = 600, width = 8, height = 6, units = "in")

smy.res2 <- smy.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 5 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps") & SampleSize %in% paste0("N = ", c(500, 1000, 2500, 5000)))

ggplot(smy.res2, aes(DIRECTION, mean_OR, fill = SampleSize)) +
    geom_col(position = "dodge", color = "black") +
    geom_hline(yintercept = med_comp, color = "black", size = 1, linetype = "dashed") +
    geom_errorbar(aes(ymin = -2 * sd_OR + mean_OR, ymax = 2 * sd_OR + mean_OR), position = position_dodge2(width = .5, padding = 0.5)) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    scale_color_grey() +
    labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI",
         x = "Intended Direction of Bias due to Missing Data",
         y = "Odds Ratio")

ggsave2(paste0(plot_filepath, "odds_ratio_column.png"), dpi = 600, width = 8, height = 6, units = "in")

# Plot bias
ggplot(small_smpl.res, aes(DIRECTION, OR_BIAS, fill = SampleSize)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1, outlier.shape = 1) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Stat. Bias of the Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Statistical Bias")

ggsave2(paste0(plot_filepath, "or_bias.png"), dpi = 600, width = 8, height = 6, units = "in")

# Plot missing cases proportions

# Filter out Complete Data Analyses
sim.res.nc <- small_smpl.res %>% filter(ANALYSIS != "COMP")

ggplot(sim.res.nc, aes(DIRECTION, PROP_MISS, fill = SampleSize)) +
    geom_boxplot(alpha = 0.75) +
    facet_wrap(MISS_TYPE~.) +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Proportion of Incomplete Cases", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Proportion of Incomplete Cases")

ggsave2(paste0(plot_filepath, "prop_miss.png"), dpi = 600, width = 8, height = 6, units = "in")

ggplot(sim.res.nc, aes(DIRECTION, OR_DIFF, fill = SampleSize)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1, outlier.shape = 1) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_y") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Difference in Race Effect Estimates from CCA/MI to the Complete Data Estimate", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Difference in the OR")

ggsave2(paste0(plot_filepath, "or_diff.png"), dpi = 600, width = 8, height = 6, units = "in")

# # Large Sample Size Plots

plot_filepath <- "Plots/sim_res_large_sample_sizes_"
large_smpl.res <- sim.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 5 imps", "MI - 8 imps") & SampleSize %in% paste0("N = ", c(10000, 25000, 50000)))

# Odds Ratio
ggplot(large_smpl.res, aes(DIRECTION, ODDS_RATIO, fill = SampleSize)) +
    geom_hline(yintercept = med_comp, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1, outlier.shape = 1) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_x", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Odds Ratio")

ggsave2(paste0(plot_filepath, "odds_ratio.png"), dpi = 600, width = 8, height = 6, units = "in")

smy.res2 <- smy.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 5 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps") & SampleSize %in% paste0("N = ", c(10000, 25000, 50000)))

ggplot(smy.res2, aes(DIRECTION, mean_OR, fill = SampleSize)) +
    geom_col(position = "dodge", color = "black") +
    geom_hline(yintercept = med_comp, color = "black", size = 1, linetype = "dashed") +
    geom_errorbar(aes(ymin = -2 * sd_OR + mean_OR, ymax = 2 * sd_OR + mean_OR), position = position_dodge2(width = .5, padding = 0.5)) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    scale_color_grey() +
    labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI",
     x = "Intended Direction of Bias due to Missing Data",
     y = "Odds Ratio")

ggsave2(paste0(plot_filepath, "odds_ratio_column.png"), dpi = 600, width = 8, height = 6, units = "in")

# Plot bias
ggplot(large_smpl.res, aes(DIRECTION, OR_BIAS, fill = SampleSize)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1, outlier.shape = 1) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Stat. Bias of the Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Statistical Bias")

ggsave2(paste0(plot_filepath, "or_bias.png"), dpi = 600, width = 8, height = 6, units = "in")

# Plot missing cases proportions

# Filter out Complete Data Analyses
sim.res.nc <- large_smpl.res %>% filter(ANALYSIS != "COMP")

ggplot(sim.res.nc, aes(DIRECTION, PROP_MISS, fill = SampleSize)) +
    geom_boxplot(alpha = 0.75) +
    facet_wrap(MISS_TYPE~.) +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Proportion of Incomplete Cases", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Proportion of Incomplete Cases")

ggsave2(paste0(plot_filepath, "prop_miss.png"), dpi = 600, width = 8, height = 6, units = "in")

ggplot(sim.res.nc, aes(DIRECTION, OR_DIFF, fill = SampleSize)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1, outlier.shape = 1) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_y") +
    theme_bw() +
    scale_fill_grey() +
    labs(title = "Dist. of the Difference in Race Effect Estimates from CCA/MI to the Complete Data Estimate", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Difference in the OR")

ggsave2(paste0(plot_filepath, "or_diff.png"), dpi = 600, width = 8, height = 6, units = "in")


# All Sample Size Plots
plot_filepath <- "Plots/sim_res_all_sample_sizes_"
all_smpl.res <- sim.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps") & !(SampleSize == "N = 500" | SampleSize == "N = 1000"))

# Odds Ratio
ggplot(all_smpl.res, aes(DIRECTION, ODDS_RATIO, color = SampleSize)) +
    geom_hline(yintercept = med_comp, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1, outlier.shape = 1) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
    theme_bw() +
    theme(legend.position = "bottom") +
    # scale_fill_grey() +
    ggthemes::scale_color_colorblind() +
    labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Odds Ratio")

all_smpl.res.mar <- sim.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps") & !(SampleSize == "N = 500" | SampleSize == "N = 1000") & MISS_TYPE == "MAR")

# Odds Ratio
ggplot(all_smpl.res.mar, aes(DIRECTION, ODDS_RATIO, color = SampleSize)) +
    geom_hline(yintercept = med_comp, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1) +
    # facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
    facet_wrap(ANALYSIS2~., scales = "free") +
    theme_bw() +
    theme(legend.position = "bottom") +
    # scale_fill_grey() +
    ggthemes::scale_color_colorblind() +
    labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Odds Ratio")


ggsave2(paste0(plot_filepath, "odds_ratio.png"), dpi = 600, width = 8, height = 6, units = "in")

smy.res2 <- smy.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 5 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps"))
ggplot(smy.res2, aes(DIRECTION, mean_OR, fill = SampleSize)) +
    geom_col(position = "dodge", color = "black") +
    geom_hline(yintercept = med_comp, color = "black", size = 1, linetype = "dashed") +
    geom_errorbar(aes(ymin = -2 * sd_OR + mean_OR, ymax = 2 * sd_OR + mean_OR), position = position_dodge2(width = .5, padding = 0.5)) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
    theme_bw() +
    # scale_fill_grey() +
    ggthemes::scale_fill_colorblind() +
    # scale_color_grey() +
    ggthemes::scale_color_colorblind() +
    labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI",
         x = "Intended Direction of Bias due to Missing Data",
         y = "Odds Ratio")

ggsave2(paste0(plot_filepath, "odds_ratio_column.png"), dpi = 600, width = 8, height = 6, units = "in")

# Plot bias
ggplot(all_smpl.res, aes(DIRECTION, OR_BIAS, color = SampleSize)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
    theme_bw() +
    # scale_fill_grey() +
    ggthemes::scale_color_colorblind() +
    labs(title = "Dist. of the Stat. Bias of the Race Effect Estimates with Complete Data, CCA, and MI", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Statistical Bias")

ggsave2(paste0(plot_filepath, "or_bias.png"), dpi = 600, width = 8, height = 6, units = "in")

# Plot missing cases proportions

# Filter out Complete Data Analyses
sim.res.nc <- all_smpl.res %>% filter(ANALYSIS != "COMP")

ggplot(sim.res.nc, aes(DIRECTION, PROP_MISS, color = SampleSize)) +
    geom_boxplot(alpha = 0.75) +
    facet_wrap(MISS_TYPE~.) +
    theme_bw() +
    # scale_fill_grey() +
    ggthemes::scale_color_colorblind() +
    labs(title = "Dist. of the Proportion of Incomplete Cases", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Proportion of Incomplete Cases")

ggsave2(paste0(plot_filepath, "prop_miss.png"), dpi = 600, width = 8, height = 6, units = "in")

ggplot(sim.res.nc, aes(DIRECTION, OR_DIFF, color = SampleSize)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(alpha = .75, outlier.shape = 1) +
    facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free_y") +
    theme_bw() +
    # scale_fill_grey() +
    ggthemes::scale_color_colorblind() +
    labs(title = "Dist. of the Difference in Race Effect Estimates from CCA/MI to the Complete Data Estimate", 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Difference in the OR")

ggsave2(paste0(plot_filepath, "or_diff.png"), dpi = 600, width = 8, height = 6, units = "in")

