
result_plots <- function(sim.res, Q = 225, N = 1000, horizontal = T) {
    # options(dplyr.summarise.inform = FALSE)
    # tbls <- result_tables(sim.res)
    # mean_OR <- tbls$ODDS_RATIO
    # mean_PM <- tbls$PROP_MISS
    # mean_DIFF <- tbls$DIFF
    # mean_BIAS <- tbls$BIAS
    
    or_plt <- my_gg_box(sim.res, X = "ODDS_RATIO", 
                        main = "Dist. of Odds Ratios of the Race Effect",
                        xlab = "Odds Ratio", Q = Q, N = N, horizontal = horizontal)
    
    bias_plt <- my_gg_box(sim.res, X = "OR_BIAS",
                          main = "Densities of Bias in Est. OR of the Race Effect", 
                          xlab = TeX("$OR_{miss} - e^\\beta$"), Q = Q, N = N, horizontal = horizontal)
    
    # Removing the complete data estimates from the table since they aren't needed for the other plots
    sim.res <- sim.res[which(sim.res$ANALYSIS != "COMP"),]
    
    pm_plt <- my_gg_box(sim.res, X = "PROP_MISS", 
                        main = "Dist. of Prop. of Missing Cases",
                        xlab = "Missing Proportion", Q = Q, N = N, horizontal = horizontal)
    
    diff_plt <- my_gg_box(sim.res, X = "OR_DIFF", 
                          main = "Densities of Diff in Est. OR of the Race Effect",
                          xlab = TeX("$OR_{miss} - OR_{comp}$"), Q = Q, N = N, horizontal = horizontal)
    
    return(list("OR" = or_plt, "PM" = pm_plt,
                "DIFF" = diff_plt, "BIAS" = bias_plt))
}

present_figs <- function(N = 500, Q = 225, m = c(3,5,8), horizontal = horizontal, ...) {
    
    sim.res <- load_sim_results(Q = Q, N = N, m = m)
    
    plt.list <- result_plots(sim.res, Q = Q, N = N, horizontal = horizontal)
    
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
        ggthemes::scale_fill_colorblind() +
        ggthemes::scale_color_colorblind() + 
        labs(y = "Percent")
    return(p1)
}

my_gg_box <- function(dat, X, main = "", xlab = "", Q = 225, N = 500, horizontal = T, ...) {
    p1 <- ggplot(dat, aes_string(x = X)) +
        geom_vline(xintercept = 0, color = "darkgray") +
        geom_boxplot(aes(y = IMPUTATION, fill = DIRECTION)) +
        # geom_violin(aes(y = IMPUTATION, fill = DIRECTION),
        #             draw_quantiles = c(0.25, 0.5, 0.75)) +
        labs(title = main,
             subtitle = paste0(Q, " iterations, n = ", N),
             x = xlab) +
        facet_grid(MISS_TYPE ~ ANALYSIS, scales = "free_y") +
        ggthemes::scale_fill_colorblind() +
        ggthemes::scale_color_colorblind() 
    
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
#              subtitle = paste0(Q, " iterations, n = ", N, ", N = ", m),
#              x = TeX("$OR_{miss} - OR_{comp}$")) +
#         facet_grid(MISS_TYPE ~ ., scales = "free_y") +
#         scale_color_colorblind()
#     return(p1)
# }