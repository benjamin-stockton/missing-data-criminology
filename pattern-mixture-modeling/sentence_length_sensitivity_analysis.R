library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(ggplot2, warn.conflicts = FALSE, quietly = TRUE)
library(mice, warn.conflicts = FALSE, quietly = TRUE)
library(GLMMadaptive)
library(broom.mixed)
library(brms)
library(VIM)
library(bayesplot)
library(loo)
setwd("pattern-mixture-modeling")

theme_set(theme_bw())
color_scheme_set("brightblue")

df <- readr::read_csv("../Data/PCS-most-serious-sentence-2010-2019-pmm.csv",
                      show_col_types = FALSE)

county_ogs <- df |>
    group_by(COUNTY) |>
    summarize(
        COUNTY_OGS = mean(OGS, na.rm = TRUE),
        NCASES = n()
    )
county_ogs

df <- df |> mutate(
    YEAR = as.factor(YEAR),
    INCAR = case_when(
        JP_MIN == 0 ~ 0,
        JP_MIN > 0 ~ 1,
        TRUE ~ NA
    ),
    SEX = case_when(
        MALE == 1 ~ "Male",
        MALE == 0 ~ "Female",
        TRUE ~ NA
    ),
    OFF_RACE = case_when(
        OFF_RACER == 1 ~ "WHITE",
        OFF_RACER == 2 ~ "BLACK",
        OFF_RACER == 3 ~ "LATINO",
        OFF_RACER == 4 ~ "OTHER",
        TRUE ~ NA
    ),
    PRVREC = case_when(
        PRSR == 0 ~ "0",
        PRSR == 1 ~ "1/2/3",
        PRSR == 2 ~ "4/5",
        PRSR == 3 ~ "REVOC/RFEL",
        TRUE ~ NA
    ),
    CRIME = case_when(
        CRIMETYPE == 1 ~ "Persons",
        CRIMETYPE == 2 ~ "Property",
        CRIMETYPE == 3 ~ "Drug",
        CRIMETYPE == 4 ~ "DUI",
        CRIMETYPE == 5 ~ "Other",
        TRUE ~ NA
    ),
    COUNTYTYPE = case_when(
        COUNTY %in% c("Allegheny", "Philadelphia") ~ "Urban",
        COUNTY %in% c("Chester", "Montgomery", "Berks", "Dauphin", "Bucks", "Lancaster", "York", "Delaware", "Northampton", "Luzerne", "Lackawanna", "Westmoreland", "Lehigh", "Erie") ~ "Medium",
        TRUE ~ "Rural"
    ),
    JP_MIN_NZ = case_when(
        JP_MIN == 0 ~ NA,
        is.na(JP_MIN) ~ NA,
        TRUE ~ JP_MIN
    )
) |>
    select(
        -c(OFF_RACER, PRSR, MALE, CRIMETYPE)
    ) |>
    left_join(county_ogs, by = "COUNTY")

df$OFF_RACE <- factor(df$OFF_RACE, levels = c("WHITE", "BLACK", "LATINO", "OTHER"))

df |> 
    select(
        JP_MIN, OFF_RACE, DOSAGE, RECMIN, PRVREC, CRIME, TRIAL, SEX
    ) |> 
    aggr(sortby = "JP_MIN", plot = FALSE) |>
    plot(numbers = TRUE, prop = FALSE)

df[,c("OGS", "DOSAGE")] <- scale(df[,c("OGS", "DOSAGE")], center = TRUE, scale = TRUE)

df[,"DOSAGEQ"] <- df[,"DOSAGE"]^2
df[,"OGSQ"] <- df[,"OGS"]^2

s <- sample(1:nrow(df), size = 10000, replace = FALSE)

df[s,] |> 
    mutate(
        jp_min_l1p = log1p(JP_MIN),
        OGSCRIMEDUI = case_when(CRIME == "DUI" ~ OGS,
                                TRUE ~ 0),
        OGSCRIMEOther = case_when(CRIME == "Other" ~ OGS,
                                  TRUE ~ 0),
        OGSCRIMEPersons = case_when(CRIME == "Persons" ~ OGS,
                                    TRUE ~ 0),
        OGSCRIMEProperty = case_when(CRIME == "Property" ~ OGS,
                                     TRUE ~ 0),
        OGSPREV123 = case_when(PRVREC == "1/2/3" ~ OGS,
                               TRUE ~ 0),
        OGSPREV45 = case_when(PRVREC == "4/5" ~ OGS,
                              TRUE ~ 0),
        OGSPREVRFEL = case_when(PRVREC == "REVOC/RFEL" ~ OGS,
                                TRUE ~ 0),
        
    ) |>
    select(starts_with("OGS"), jp_min_l1p) |>
    GGally::ggpairs()

df |>
    group_by(CRIME, PRVREC) |>
    summarize(
        mean_sen_len = mean(JP_MIN, na.rm = TRUE),
        sd_sen_len = sd(JP_MIN, na.rm = TRUE),
        p_incar = mean(JP_MIN == 0, na.rm = TRUE)
    ) |>
    print(n = 21)

imps0 <- mice(df[s,], m = 1, method = "pmm", maxit = 0)

mthd <- imps0$method
mthd["JP_MIN_MON"] <-  "~I(JP_MIN / 30)"
mthd["INCAR"] <- "~I(ifelse(JP_MIN > 0, 1, 0))"
# mthd[c("RECMIN", "INCAR", "OFF_RACE", "PRVREC")] <- "cart"
mthd

pred_mat <- imps0$predictorMatrix
pred_mat[,"JPR_ID"] <- 0
pred_mat

imps <- mice(df, m = 5, method = mthd, predictorMatrix = pred_mat, maxit = 10)

plot(imps)

saveRDS(imps, "fits/mice-rf-subsample.rds")

imps <- readRDS("fits/mice-rf-subsample.rds")
imp_list <- complete(imps, "all", include = FALSE)


# Write function to obtain estimate under given delta parameter
scale_MNAR_f <- function(scale){

    df2 <- df |> 
        select(
            JP_MIN_NZ, DOSAGE, DOSAGEQ, OGS, OGSQ, RECMIN, TRIAL,
            YEAR, COUNTY, INCAR, SEX, OFF_RACE, PRVREC, CRIME, 
            COUNTYTYPE, COUNTY_OGS, NCASES
        )

    imps0 <- mice(df2[s,], m = 1, method = "pmm", maxit = 0)
    
    mthd <- imps0$method
    # mthd["JP_MIN_MON"] <-  "~I(JP_MIN / 30)"
    mthd["INCAR"] <- "~I(ifelse(JP_MIN_NZ > 0, 1, 0))"
    # mthd[c("RECMIN", "INCAR", "OFF_RACE", "PRVREC")] <- "cart"
    mthd
    
    pred_mat <- imps0$predictorMatrix
    # pred_mat[,"JPR_ID"] <- 0
    pred_mat
    
    MAR_imps <- mice(df2, m = 5, method = mthd, predictorMatrix = pred_mat, maxit = 5)
    
    MNAR_long <- MAR_imps |>
        complete(action='long', include=TRUE) |>
        mutate(JP_MIN2 = case_when(
            INCAR == 0 ~ JP_MIN_NZ * scale,
            TRUE ~ JP_MIN_NZ),
            log_JP_MIN2 = log(JP_MIN2))
    
    MNAR_imps <- as.mids(MNAR_long)
    
    fit_scale <- with(MNAR_imps, lme4::lmer(log_JP_MIN2 ~ DOSAGE + DOSAGEQ + SEX + OFF_RACE + PRVREC * OGS + OGSQ + PRVREC * CRIME + RECMIN + TRIAL + (1 + COUNTY_OGS + NCASES || COUNTY)))
    
    
    params_scale <- summary(pool(fit_scale)) |> 
        filter(term == "OFF_RACEBLACK") 
    
    desc_scale <- with(MNAR_imps, expr=c("Y_mean"=mean(log_JP_MIN2), "Y_sd"=stats::sd(log_JP_MIN2)))
    Q_mean <- unlist(lapply(desc_scale$analyses, function(vec) {return(vec[1])}))
    U_mean <- unlist(lapply(desc_scale$analyses, function(vec) {return(vec[2]^2)}))
    desc_scale <- pool.scalar(Q_mean, U_mean, n = nrow(df)) |>
        as.data.frame() |>
        select(qhat, t)
    colnames(desc_scale) <- c("Y_mean", "Y_var")
    
    cbind(params_scale, desc_scale[1,])
    
}

# Set range of plausible scale values
sens_scale <- c(0.5, 1, 1.5)

# Map output
output_scale <- purrr::map_df(sens_scale, scale_MNAR_f)

rownames(output_scale) <- c()

output_scale <- cbind(sens_scale, output_scale)

output_scale |>
    select(-term) |>
    mutate(across(c("estimate":"std.error"), round, 3),
           across(c("statistic":"df", "Y_mean":"Y_var"), round, 2))

output_scale |>
    ggplot(aes(x = sens_scale, y = estimate)) +
    geom_smooth() + 
    geom_point() + 
    geom_hline(yintercept=0, linetype="dashed", color = "red") +
    theme_minimal()
