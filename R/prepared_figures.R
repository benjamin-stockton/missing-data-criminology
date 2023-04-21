source("R/helpers.R")
source("R/results_figures.R")
source("R/results_tables.R")
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
    tmp <- load_sim_results(Q = Q, N = n, m = M)

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

#####################################################
# Plotting Results
#####################################################

prep_figs <- function(data, plot_filepath, fig.width = 8, fig.height = 6) {
    p1 <- ggplot(data) +
        # ggh4x::facet_grid2(vars(MISS_TYPE), vars(ANALYSIS2), 
                           # scales = "free", space = "free", independent = "y") +
        facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
        theme_bw() +
        ggthemes::scale_fill_colorblind() +
        ggthemes::scale_color_colorblind() +
        theme(legend.position = "bottom")
    
    # Odds Ratio
    or <- p1 + 
        geom_hline(yintercept = med_comp, color = "gray") +
        geom_boxplot(aes(DIRECTION, ODDS_RATIO, color = SampleSize), alpha = .75, outlier.shape = 1) +
        labs(title = "Dist. of Race Effect Estimates with Complete Data, CCA, and MI", 
             x = "Intended Direction of Bias due to Missing Data",
             y = "Odds Ratio")
    
    ggsave2(paste0(plot_filepath, "odds_ratio.tiff"), dpi = 600, width = 8, height = 6, units = "in")
    
    # Plot bias
    bias <- p1 +
        geom_hline(yintercept = 0, color = "gray") +
        geom_boxplot(aes(DIRECTION, OR_BIAS, color = SampleSize), alpha = .75, outlier.shape = 1) +
        labs(title = "Dist. of the Stat. Bias of the Race Effect Estimates with Complete Data, CCA, and MI", 
             x = "Intended Direction of Bias due to Missing Data",
             y = "Statistical Bias")
    
    ggsave2(paste0(plot_filepath, "or_bias.tiff"), dpi = 600, width = 8, height = 6, units = "in")
    
    # Plot missing cases proportions
    
    # Filter out Complete Data Analyses
    sim.res.nc <- data %>% filter(ANALYSIS != "COMP")
    p2 <- ggplot(sim.res.nc) +
        # ggh4x::facet_grid2(vars(MISS_TYPE), vars(ANALYSIS2),
        # scales = "free", space = "free", independent = "y") +
        facet_grid(MISS_TYPE ~ ANALYSIS2, scales = "free", space = "free") +
        theme_bw() +
        ggthemes::scale_fill_colorblind() +
        ggthemes::scale_color_colorblind() +
        theme(legend.position = "bottom")
    
    pmiss <- p2 + 
        geom_boxplot(aes(DIRECTION, PROP_MISS, color = SampleSize), alpha = 0.75) +
        facet_grid(MISS_TYPE~., scales = "free") +
        labs(title = "Dist. of the Proportion of Incomplete Cases", 
             x = "Intended Direction of Bias due to Missing Data",
             y = "Proportion of Incomplete Cases")
    
    ggsave2(paste0(plot_filepath, "prop_miss.tiff"), dpi = 600, width = 8, height = 6, units = "in")
    
    ordiff <- p2 + 
        geom_hline(yintercept = 0, color = "gray") +
        geom_boxplot(aes(DIRECTION, OR_DIFF, color = SampleSize), alpha = .75, outlier.shape = 1) +
        labs(title = "Dist. of the Difference in Race Effect Estimates from CCA/MI to the Complete Data Estimate", 
             x = "Intended Direction of Bias due to Missing Data",
             y = "Difference in the OR")
    
    # cowplot::plot_grid(or, bias, pmiss, ordiff, nrow = 2)
    ggsave2(paste0(plot_filepath, "or_diff.tiff"), dpi = 600, width = 8, height = 6, units = "in")
}

# Small Sample Size Plots
ss_fp <- "Plots/sim_res_small_sample_sizes_"
small_smpl.res <- sim.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps") & SampleSize %in% paste0("N = ", c(500, 1000, 2500, 5000)))

prep_figs(data = small_smpl.res, plot_filepath = ss_fp)


# Large Sample Size Plots

ls_fp <- "Plots/sim_res_large_sample_sizes_"
large_smpl.res <- sim.res %>% filter(SampleSize %in% paste0("N = ", c(10000, 25000, 50000)))

prep_figs(data = large_smpl.res, plot_filepath = ls_fp)


# All Sample Size Plots
all_fp <- "Plots/sim_res_all_sample_sizes_"
all_smpl.res <- sim.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps"))

prep_figs(data = all_smpl.res, plot_filepath = all_fp)

all_fp_mar <- "Plots/sim_res_all_sample_sizes_MAR_"
all_smpl.res.mar <- sim.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 5 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps") & MISS_TYPE == "MAR")

prep_figs(data = all_smpl.res.mar, plot_filepath = all_fp_mar)


all_fp_mnar <- "Plots/sim_res_all_sample_sizes_MNAR_"
all_smpl.res.mnar <- sim.res %>% filter(ANALYSIS2 %in% c("COMP", "CCA", "MI - 3 imps", "MI - 5 imps", "MI - 8 imps", "MI - 25 imps", "MI - 100 imps") & MISS_TYPE == "MNAR")

prep_figs(data = all_smpl.res.mnar, plot_filepath = all_fp_mnar)

mar.nc <- all_smpl.res.mar %>% filter(ANALYSIS2 != "COMP")
ggplot(mar.nc) +
    # ggh4x::facet_grid2(vars(MISS_TYPE), vars(ANALYSIS2),
    # scales = "free", space = "free", independent = "y") +
    # facet_grid(. ~ ANALYSIS2, scales = "free", space = "free") +
    ggh4x::facet_wrap2(vars(ANALYSIS2), nrow = 2, trim_blank = T, axes = "all", scales = "free") +
    theme_bw() +
    ggthemes::scale_fill_colorblind() +
    ggthemes::scale_color_colorblind() +
    theme(legend.position = "bottom") + 
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(aes(DIRECTION, OR_DIFF, color = SampleSize), alpha = .75, outlier.shape = 1) +
    labs(title = latex2exp::TeX("Dist. of $d_{Dir}^{Analysis} = \\hat{OR}^{Analysis}_{Race}(D_{Dir}) - \\hat{OR}_{Race}(D_{Comp})$ under MAR simulations"), 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Difference in the OR")

mnar.nc <- all_smpl.res.mnar %>% filter(ANALYSIS2 != "COMP")
ggplot(mnar.nc) +
    # ggh4x::facet_grid2(vars(MISS_TYPE), vars(ANALYSIS2),
    # scales = "free", space = "free", independent = "y") +
    # facet_grid(. ~ ANALYSIS2, scales = "free", space = "free") +
    ggh4x::facet_wrap2(vars(ANALYSIS2), nrow = 2, trim_blank = T, axes = "all", scales = "free") +
    theme_bw() +
    ggthemes::scale_fill_colorblind() +
    ggthemes::scale_color_colorblind() +
    theme(legend.position = "bottom") + 
    geom_hline(yintercept = 0, color = "gray") +
    geom_boxplot(aes(DIRECTION, OR_DIFF, color = SampleSize), alpha = .75, outlier.shape = 1) +
    labs(title = latex2exp::TeX("Dist. of $d_{Dir}^{Analysis} = \\hat{OR}^{Analysis}_{Race}(D_{Dir}) - \\hat{OR}_{Race}(D_{Comp})$ under MNAR simulations"), 
         x = "Intended Direction of Bias due to Missing Data",
         y = "Difference in the OR")

# se.beta <- all_smpl.res %>%
#     mutate(
#         rmi = case_when(MISS_TYPE == "MNAR" ~ 0.097,
#                         MISS_TYPE == "MAR" ~ 0.05)) %>%
#     group_by(MISS_TYPE, SampleSize, ANALYSIS2) %>%
#     summarize(
#         mean_rmi = mean(rmi),
#         se_beta = sd(log(ODDS_RATIO)),
#         cv_beta = 0.01 / se_beta,
#         M = 1 + 0.5 * (mean_rmi / cv_beta)^2
#     )
# se.beta %>% print(n = 78)



