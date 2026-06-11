#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: All DSEI (indigenous health subdistrict) analyses.
#          Produces Tables S2, S5, S7, S10, S11, S18, S21, S23.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/01_setup.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data and create derived variables -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dsei_df <- read_csv(paste0(data_dir, "dsei_panel_goldneighbor_0319.csv"),
                    show_col_types = FALSE)
dsei_df$seqid <- as.character(dsei_df$seqid)

# Derived outcome variables
dsei_df$falciparum_api <- 1000 * dsei_df$malaria_n_falciparum_mix / dsei_df$population_polobase
dsei_df$non_f_api      <- 1000 * dsei_df$malaria_n_not_f / dsei_df$population_polobase

cat("DSEI panel loaded:", nrow(dsei_df), "rows,", ncol(dsei_df), "cols\n")
cat("N polo bases:", length(unique(dsei_df$seqid)), "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Instrument selection: manual backward elimination from rlasso candidates -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- DSEI main (unweighted) ---
# rlasso selection (Jan 25, 2026): mp3_delta_cs pp3_alfa_bj pp3_delta_in
#   pp3_gamma_mlp c2i pp4_gamma_po (6 instruments)
cat("\n=== Instrument selection: DSEI main (unweighted) ===\n")

# Test all 6 rlasso-selected instruments
fmla_all6 <- as.formula("malaria_api_by_polobase ~ 1 | seqid + year |
    goldmine_area ~ mp3_delta_cs_2yrgp + pp3_delta_in_2yrgp +
    pp3_gamma_mlp_2yrgp + c2i_2yrgp + pp4_gamma_po_2yrgp + pp3_alfa_bj_2yrgp")
iv_all6 <- feols(fmla_all6, data = dsei_df, cluster = "seqid")
cat("All 6 IVs — first stage:\n")
print(summary(iv_all6, stage = 1)$coeftable)

# Test subset for collinearity: drop pp3_alfa_bj and pp4_gamma_po
fmla_4iv <- as.formula("malaria_api_by_polobase ~ 1 | seqid + year |
    goldmine_area ~ mp3_delta_cs_2yrgp + pp3_delta_in_2yrgp +
    pp3_gamma_mlp_2yrgp + c2i_2yrgp")
iv_4iv <- feols(fmla_4iv, data = dsei_df, cluster = "seqid")
cat("\n4 IVs (drop collinear) — first stage:\n")
print(summary(iv_4iv, stage = 1)$coeftable)
cat("F-stat (4 IVs):", round(fitstat(iv_4iv, type = "ivf")[[1]]$stat, 1), "\n")

# Remove least significant: keep pp3_delta_in + pp3_gamma_mlp (highest F-stat)
fmla_2iv <- as.formula("malaria_api_by_polobase ~ 1 | seqid + year |
    goldmine_area ~ pp3_delta_in_2yrgp + pp3_gamma_mlp_2yrgp")
iv_2iv <- feols(fmla_2iv, data = dsei_df, cluster = "seqid")
cat("\nFinal 2 IVs — first stage:\n")
print(summary(iv_2iv, stage = 1)$coeftable)
cat("F-stat (2 IVs):", round(fitstat(iv_2iv, type = "ivf")[[1]]$stat, 1), "\n")
cat("Final selection: pp3_delta_in_2yrgp + pp3_gamma_mlp_2yrgp\n")

# --- DSEI weighted (Table S18) ---
# rlasso weighted selection: mp3_delta_cs pp3_alfa_bj pp3_delta_in
#   pp3_gamma_mlp pp3cc c2i pp4_gamma_po (7 instruments)
cat("\n=== Instrument selection: DSEI weighted ===\n")

# Test all 7 rlasso-selected instruments (weighted)
fmla_w_all7 <- as.formula("malaria_api_by_polobase ~ 1 | seqid + year |
    goldmine_area ~ mp3_delta_cs_2yrgp + pp3_delta_in_2yrgp +
    pp3_gamma_mlp_2yrgp + c2i_2yrgp + pp4_gamma_po_2yrgp +
    pp3_alfa_bj_2yrgp + pp3cc_2yrgp")
iv_w_all7 <- feols(fmla_w_all7, data = dsei_df,
                   weights = ~population_polobase, cluster = "seqid")
cat("All 7 IVs (weighted) — first stage:\n")
print(summary(iv_w_all7, stage = 1)$coeftable)

# Drop collinear vars (pp3_gamma_mlp, pp4_gamma_po)
fmla_w_5iv <- as.formula("malaria_api_by_polobase ~ 1 | seqid + year |
    goldmine_area ~ mp3_delta_cs_2yrgp + pp3_delta_in_2yrgp +
    c2i_2yrgp + pp3_alfa_bj_2yrgp + pp3cc_2yrgp")
iv_w_5iv <- feols(fmla_w_5iv, data = dsei_df,
                  weights = ~population_polobase, cluster = "seqid")
cat("\n5 IVs (drop collinear) — first stage:\n")
print(summary(iv_w_5iv, stage = 1)$coeftable)

# Remove most insignificant (pp3cc)
fmla_w_4iv <- as.formula("malaria_api_by_polobase ~ 1 | seqid + year |
    goldmine_area ~ mp3_delta_cs_2yrgp + pp3_delta_in_2yrgp +
    c2i_2yrgp + pp3_alfa_bj_2yrgp")
iv_w_4iv <- feols(fmla_w_4iv, data = dsei_df,
                  weights = ~population_polobase, cluster = "seqid")
cat("\nFinal 4 IVs (weighted) — first stage:\n")
print(summary(iv_w_4iv, stage = 1)$coeftable)
cat("F-stat (4 IVs):", round(fitstat(iv_w_4iv, type = "ivf")[[1]]$stat, 1), "\n")
cat("Final selection: mp3_delta_cs + pp3_delta_in + c2i + pp3_alfa_bj\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S2: Main IV results (indigenous population) -----
# Output: malaria_dsei_neighbor_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S2: Main IV results (DSEI) ===\n")

dsei_outcomes <- c("malaria_api_by_polobase", "falciparum_api", "non_f_api")

# Run IV: no covariates + select covariates for each outcome
models_s2 <- list()
for (y in dsei_outcomes) {
  models_s2[[length(models_s2) + 1]] <- feols(build_dsei_iv_formula(y, "none"),
                                              data = dsei_df, cluster = "seqid")
  models_s2[[length(models_s2) + 1]] <- feols(build_dsei_iv_formula(y, "select"),
                                              data = dsei_df, cluster = "seqid")
}
names(models_s2) <- c(
  "Malaria API", "Malaria API",
  "Falciparum API", "Falciparum API",
  "Non-falciparum API", "Non-falciparum API"
)

# Mean of DV row
mean_malaria_dsei    <- round(mean(dsei_df$malaria_api_by_polobase, na.rm = TRUE), 2)
mean_falciparum_dsei <- round(mean(dsei_df$falciparum_api, na.rm = TRUE), 2)
mean_non_f_dsei      <- round(mean(dsei_df$non_f_api, na.rm = TRUE), 2)

mean_row_s2 <- data.frame(
  term = "Mean of DV",
  m1 = mean_malaria_dsei, m2 = mean_malaria_dsei,
  m3 = mean_falciparum_dsei, m4 = mean_falciparum_dsei,
  m5 = mean_non_f_dsei, m6 = mean_non_f_dsei
)

out_s2 <- paste0(table_dir, "malaria_dsei_neighbor_area.tex")
modelsummary(models_s2, stars = TRUE,
             coef_rename = coef_rename_dsei,
             gof_map = gof_map_iv_dsei,
             add_rows = mean_row_s2,
             output = out_s2)
strip_table_float(out_s2)
add_spanning_headers(out_s2)
cat("Saved:", out_s2, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S7: First-stage diagnostics (DSEI) -----
# Output: first_stage_diagnostics_dsei_neighbor_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S7: First-stage diagnostics (DSEI) ===\n")

# Use the with-covariates models (indices 2, 4, 6)
dsei_cov_models <- models_s2[c(2, 4, 6)]
spec_names_s7 <- c("Malaria API", "Falciparum API", "Non-falciparum API")

out_s7 <- paste0(table_dir, "first_stage_diagnostics_dsei_neighbor_area.tex")
make_iv_summary_table(dsei_cov_models, spec_names_s7, output_path = out_s7)
strip_table_float(out_s7)
cat("Saved:", out_s7, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S10: DSEI exclusion restriction — cross-sectional rock shares -----
# Output: dsei_exclusion_rockshare_2003_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S10: Exclusion restriction (cross-sectional OLS, DSEI) ===\n")

dsei_03 <- dsei_df %>% filter(year == 2003)

excl_outcomes_dsei <- c("all_ag", "soy", "forest", "pasture",
                        "mean_temp", "annual_total_precip")
excl_names_dsei <- c("Agriculture", "Soy", "Forest", "Pasture",
                     "Avg. temp.", "Precip.")

models_s10 <- lapply(excl_outcomes_dsei, function(y) {
  lm(as.formula(paste0(y, " ~ pp3_delta_in + pp3_gamma_mlp")), data = dsei_03)
})
names(models_s10) <- excl_names_dsei

out_s10 <- paste0(table_dir, "dsei_exclusion_rockshare_2003_area.tex")
modelsummary(models_s10, stars = TRUE,
             statistic = NULL,
             coef_rename = c(
               "pp3_delta_in" = "Ingarana",
               "pp3_gamma_mlp" = "Granito Pepita"
             ),
             gof_map = gof_map_iv_dsei,
             output = out_s10)
strip_table_float(out_s10)
cat("Saved:", out_s10, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S11: DSEI exclusion restriction — IV first-stage on covariates -----
# Output: dsei_exclusion_iv_firststage.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S11: Exclusion restriction (IV first-stage on covariates, DSEI) ===\n")

excl_iv_outcomes_dsei <- c("all_ag", "soy", "forest", "pasture",
                           "mean_temp", "annual_total_precip")

models_s11 <- lapply(excl_iv_outcomes_dsei, function(cov_var) {
  fmla <- as.formula(paste0(
    "malaria_api_by_polobase ~ 1 | seqid + year | ",
    cov_var, " ~ ", dsei_instruments
  ))
  m <- feols(fmla, data = dsei_df, cluster = "seqid")
  summary(m, stage = 1)
})
names(models_s11) <- excl_names_dsei

out_s11 <- paste0(table_dir, "dsei_exclusion_iv_firststage.tex")
modelsummary(models_s11, stars = TRUE,
             statistic = NULL,
             coef_rename = c(
               "fit_goldmine_area" = "Gold mining area",
               "pp3_delta_in_2yrgp" = "Ingarana x 2 yr gold price",
               "pp3_gamma_mlp_2yrgp" = "Granito Pepita x 2 yr gold price"
             ),
             gof_map = gof_map_iv_dsei,
             output = out_s11)
strip_table_float(out_s11)
cat("Saved:", out_s11, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S5: Balance table (DSEI, 2003) -----
# Output: balance_table_dsei_2003.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S5: Balance table (DSEI, 2003) ===\n")

# Define gold-producing vs neighbor status
gm_dsei <- dsei_df %>%
  select(seqid, goldmine_area) %>%
  group_by(seqid) %>%
  summarise(total_gold = sum(goldmine_area), .groups = "drop") %>%
  mutate(dsei_status = ifelse(total_gold > 0, "gold_producing", "gold_neighbor")) %>%
  select(-total_gold)

dsei_df_status <- dsei_df %>% left_join(gm_dsei, by = "seqid")

dsei_2003 <- dsei_df_status %>% filter(year == 2003)

balance_vars_dsei <- c("population_polobase", "area_km2", "mean_temp",
                       "annual_total_precip",
                       "percent_forest", "percent_soy", "percent_all_ag",
                       "percent_pasture")

balancetbl_dsei <- dsei_2003 %>%
  mutate(dsei_status = group_labels_dsei[dsei_status]) %>%
  select(dsei_status, all_of(balance_vars_dsei)) %>%
  tbl_summary(
    by = dsei_status,
    label = as.list(var_labels_dsei[intersect(names(var_labels_dsei), balance_vars_dsei)]),
    statistic = all_continuous() ~ "{mean} ({sd})",
    missing = "no"
  ) %>%
  add_p(test = all_continuous() ~ "t.test") %>%
  modify_header(
    label ~ "**Variable**",
    all_stat_cols() ~ "**{level}** (N = {n})"
  )

out_s5 <- paste0(table_dir, "balance_table_dsei_2003.tex")
balancetbl_dsei %>%
  as_gt() %>%
  gt::gtsave(filename = out_s5)
strip_table_float(out_s5)
cat("Saved:", out_s5, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S18: Weighted IV results (DSEI) -----
# Output: malaria_dsei_neighbor_weighted_area.tex
# Uses DIFFERENT instruments from main analysis (4 IVs from weighted LASSO)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S18: Weighted IV results (DSEI) ===\n")

# dsei_weighted_instruments defined in 01_setup.R (from backward elimination above)

# Weighted covariates
dsei_weighted_select_cov <- "population_polobase + mean_temp + annual_total_precip + forest"

fmla_w_nocov <- as.formula(paste0(
  "malaria_api_by_polobase ~ 1 | seqid + year | goldmine_area ~ ",
  dsei_weighted_instruments
))
w_nocov <- feols(fmla_w_nocov, data = dsei_df,
                 weights = ~ population_polobase, cluster = "seqid")

fmla_w_cov <- as.formula(paste0(
  "malaria_api_by_polobase ~ ", dsei_weighted_select_cov,
  " | seqid + year | goldmine_area ~ ", dsei_weighted_instruments
))
w_cov <- feols(fmla_w_cov, data = dsei_df,
               weights = ~ population_polobase, cluster = "seqid")

models_s18 <- list("Malaria API" = w_nocov, "Malaria API" = w_cov)

mean_dsei_w <- round(mean(dsei_df$malaria_api_by_polobase, na.rm = TRUE), 2)
mean_row_s18 <- data.frame(term = "Mean of DV", m1 = mean_dsei_w, m2 = mean_dsei_w)

out_s18 <- paste0(table_dir, "malaria_dsei_neighbor_weighted_area.tex")
modelsummary(models_s18, stars = TRUE,
             coef_rename = coef_rename_dsei,
             gof_map = gof_map_iv_dsei,
             add_rows = mean_row_s18,
             output = out_s18)
strip_table_float(out_s18)
cat("Saved:", out_s18, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S21: DSEI results excluding 2019 -----
# Output: malaria_dsei_neighbor_area_drop2019.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S21: DSEI results excluding 2019 ===\n")

dsei_no2019 <- dsei_df %>% filter(year != 2019)
cat("Sample after dropping 2019:", nrow(dsei_no2019), "observations\n")

# Derived variables for filtered sample
dsei_no2019$falciparum_api <- 1000 * dsei_no2019$malaria_n_falciparum_mix / dsei_no2019$population_polobase
dsei_no2019$non_f_api      <- 1000 * dsei_no2019$malaria_n_not_f / dsei_no2019$population_polobase

models_s21 <- list()
for (y in dsei_outcomes) {
  models_s21[[length(models_s21) + 1]] <- feols(build_dsei_iv_formula(y, "none"),
                                                data = dsei_no2019, cluster = "seqid")
  models_s21[[length(models_s21) + 1]] <- feols(build_dsei_iv_formula(y, "select"),
                                                data = dsei_no2019, cluster = "seqid")
}
names(models_s21) <- c(
  "Malaria API", "Malaria API",
  "Falciparum API", "Falciparum API",
  "Non-falciparum API", "Non-falciparum API"
)

mean_malaria_d19    <- round(mean(dsei_no2019$malaria_api_by_polobase, na.rm = TRUE), 2)
mean_falciparum_d19 <- round(mean(dsei_no2019$falciparum_api, na.rm = TRUE), 2)
mean_non_f_d19      <- round(mean(dsei_no2019$non_f_api, na.rm = TRUE), 2)

mean_row_s21 <- data.frame(
  term = "Mean of DV",
  m1 = mean_malaria_d19, m2 = mean_malaria_d19,
  m3 = mean_falciparum_d19, m4 = mean_falciparum_d19,
  m5 = mean_non_f_d19, m6 = mean_non_f_d19
)

out_s21 <- paste0(table_dir, "malaria_dsei_neighbor_area_drop2019.tex")
modelsummary(models_s21, stars = TRUE,
             coef_rename = coef_rename_dsei,
             gof_map = gof_map_iv_dsei,
             add_rows = mean_row_s21,
             output = out_s21)
strip_table_float(out_s21)
add_spanning_headers(out_s21)
cat("Saved:", out_s21, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S23: DSEI falsification tests -----
# Output: dsei_neighbor_falsification.tex
# Uses sample-specific instruments from rlasso CSV + backward elimination
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S23: DSEI falsification ===\n")

# Load DSEI health panel and merge
hp <- read.csv(paste0(data_dir, "dsei_health_panel.csv"))
hp$polo_base_lower <- tolower(hp$polo_base)
hp$dsei_lower <- tolower(hp$dsei)

merged_df <- merge(dsei_df, hp,
                   by.x = c("year", "dsei_gesta", "polo"),
                   by.y = c("year", "dsei_lower", "polo_base_lower"),
                   all.x = TRUE)

# Subset to polo bases with sufficient births for mortality outcomes
df_births5 <- merged_df %>% filter(!is.na(birth_total_live_births) & birth_total_live_births >= 5)

cat(sprintf("Full panel:                   %d obs, %d polo bases\n",
            nrow(merged_df), n_distinct(merged_df$seqid)))
cat(sprintf("Births >= 5 (for IMR, etc.):  %d obs, %d polo bases\n",
            nrow(df_births5), n_distinct(df_births5$seqid)))

# Instruments for falsification: hardcoded from backward elimination
# (rlasso selects 6 IVs; after backward elimination, 2 remain)
iv_mortality <- c("pp3_delta_in_2yrgp", "pp3_gamma_mlp_2yrgp")

make_fals_fmla <- function(outcome, ivs) {
  iv_str <- paste(ivs, collapse = " + ")
  as.formula(paste0(outcome, " ~ 1 | seqid + year | goldmine_area ~ ", iv_str))
}

# Estimate falsification models
imr    <- feols(make_fals_fmla("child_mort_imr", iv_mortality), data = df_births5, cluster = "seqid")
mort_Q <- feols(make_fals_fmla("child_mort_congenital_Q_per1k", iv_mortality), data = df_births5, cluster = "seqid")
mort_I <- feols(make_fals_fmla("child_mort_circulatory_I_per1k", iv_mortality), data = df_births5, cluster = "seqid")
mort_E <- feols(make_fals_fmla("child_mort_endocrine_E_per1k", iv_mortality), data = df_births5, cluster = "seqid")
mort_Y <- feols(make_fals_fmla("child_mort_external_VWXY_per1k", iv_mortality), data = df_births5, cluster = "seqid")

models_s23 <- list(
  "Total IMR"       = imr,
  "Congenital IMR"  = mort_Q,
  "Circulatory IMR" = mort_I,
  "Endocrine IMR"   = mort_E,
  "External IMR"    = mort_Y
)

out_s23 <- paste0(table_dir, "dsei_neighbor_falsification.tex")
modelsummary(models_s23, stars = TRUE,
             fmt = fmt_smart,
             coef_rename = coef_rename_dsei,
             gof_map = gof_map_iv_dsei,
             output = out_s23)
strip_table_float(out_s23)

# Append table note
note_lines <- c(
  "\\vspace{2pt}",
  "\\begin{minipage}{\\linewidth}",
  "\\footnotesize",
  "\\textit{Notes:} IMR = infant mortality rate.",
  "Congenital = ICD-10 chapter Q (congenital malformations).",
  "Circulatory = ICD-10 chapter I.",
  "Endocrine = ICD-10 chapter E (endocrine, nutritional, and metabolic diseases).",
  "External = ICD-10 chapters V--Y (external causes of morbidity and mortality).",
  "Mortality outcomes (IMR, Congenital, Circulatory, Endocrine, External) are per 1,000 live births.",
  "Observations with fewer than 5 live births (for mortality and stillborn outcomes) are excluded.",
  "\\end{minipage}"
)
tex <- readLines(out_s23, warn = FALSE)
tex <- c(tex[1:(length(tex) - 1)], note_lines, tex[length(tex)])
writeLines(tex, out_s23)
cat("Saved:", out_s23, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== DSEI analysis complete ===\n")
cat("Tables generated:\n")
cat("  S2:  malaria_dsei_neighbor_area.tex\n")
cat("  S5:  balance_table_dsei_2003.tex\n")
cat("  S7:  first_stage_diagnostics_dsei_neighbor_area.tex\n")
cat("  S10: dsei_exclusion_rockshare_2003_area.tex\n")
cat("  S11: dsei_exclusion_iv_firststage.tex\n")
cat("  S18: malaria_dsei_neighbor_weighted_area.tex\n")
cat("  S21: malaria_dsei_neighbor_area_drop2019.tex\n")
cat("  S23: dsei_neighbor_falsification.tex\n")
