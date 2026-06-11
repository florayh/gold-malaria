#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: All municipality-level analyses. Produces Tables S1, S3, S4, S6, S8,
#          S9, S12, S13, S17, S19, S20, S22, S24, S25, S28.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/01_setup.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data and create derived variables -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod_df <- read_csv(paste0(data_dir, "mu_panel_goldneighbor_0319.csv"),
                   show_col_types = FALSE)

# Derived outcome variables
mod_df$falciparum_api        <- 1000 * mod_df$n_falciparum_mix / mod_df$pop
mod_df$non_f_api             <- 1000 * mod_df$n_not_f / mod_df$pop
mod_df$n_malaria_death_per1000 <- 1000 * mod_df$n_malaria_deaths / mod_df$pop

cat("Municipality panel loaded:", nrow(mod_df), "rows,", ncol(mod_df), "cols\n")
cat("N municipalities:", length(unique(mod_df$cd_mun)), "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Instrument selection: manual backward elimination from rlasso candidates -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# rlasso selection (Jan 25, 2026): pp3cc pp3_delta_inc pp3_gamma_mlp (3 instruments)
# pp3cc is NOT individually significant (p=0.77) and weakens first stage
# (F-stat drops from 1002 to 668 when included)
cat("\n=== Instrument selection: municipality main ===\n")

# Test all 3 rlasso-selected instruments
fmla_fs_3iv <- as.formula("malaria_allpop_api ~ 1 | cd_mun + year |
    gold_mining_area ~ pp3cc_2yrgp + pp3_delta_inc_2yrgp + pp3_gamma_mlp_2yrgp")
fs_3iv <- feols(fmla_fs_3iv, data = mod_df, cluster = "cd_mun")
cat("All 3 IVs — first stage:\n")
print(summary(fs_3iv, stage = 1)$coeftable)
cat("F-stat (3 IVs):", round(fitstat(fs_3iv, type = "ivf")[[1]]$stat, 1), "\n")

# Drop pp3cc (insignificant, p=0.77): final 2 IVs
fmla_fs_2iv <- as.formula("malaria_allpop_api ~ 1 | cd_mun + year |
    gold_mining_area ~ pp3_delta_inc_2yrgp + pp3_gamma_mlp_2yrgp")
fs_2iv <- feols(fmla_fs_2iv, data = mod_df, cluster = "cd_mun")
cat("\nFinal 2 IVs — first stage:\n")
print(summary(fs_2iv, stage = 1)$coeftable)
cat("F-stat (2 IVs):", round(fitstat(fs_2iv, type = "ivf")[[1]]$stat, 1), "\n")
cat("Final selection: pp3_delta_inc_2yrgp + pp3_gamma_mlp_2yrgp\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S1: Main IV results (general population) -----
# Output: malaria_general_neighbormun_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S1: Main IV results ===\n")

outcomes_main <- c("malaria_allpop_api", "falciparum_api", "non_f_api",
                   "n_malaria_death_per1000")

# Run IV: no covariates + select covariates for each outcome
models_s1 <- list()
for (y in outcomes_main) {
  models_s1[[length(models_s1) + 1]] <- feols(build_mu_iv_formula(y, "none"),
                                              data = mod_df, cluster = "cd_mun")
  models_s1[[length(models_s1) + 1]] <- feols(build_mu_iv_formula(y, "select"),
                                              data = mod_df, cluster = "cd_mun")
}
names(models_s1) <- c(
  "Malaria API", "Malaria API",
  "Falciparum API", "Falciparum API",
  "Non-falciparum API", "Non-falciparum API",
  "N deaths from malaria", "N deaths from malaria"
)

# Mean of DV row
mean_malaria    <- round(mean(mod_df$malaria_allpop_api, na.rm = TRUE), 2)
mean_falciparum <- round(mean(mod_df$falciparum_api, na.rm = TRUE), 2)
mean_non_f      <- round(mean(mod_df$non_f_api, na.rm = TRUE), 2)
mean_death      <- round(mean(mod_df$n_malaria_death_per1000, na.rm = TRUE), 4)

mean_row_s1 <- data.frame(
  term = "Mean of DV",
  m1 = mean_malaria, m2 = mean_malaria,
  m3 = mean_falciparum, m4 = mean_falciparum,
  m5 = mean_non_f, m6 = mean_non_f,
  m7 = mean_death, m8 = mean_death
)

out_s1 <- paste0(table_dir, "malaria_general_neighbormun_area.tex")
modelsummary(models_s1, stars = TRUE,
             fmt = fmt_smart,
             coef_rename = coef_rename_mu,
             gof_map = gof_map_iv,
             add_rows = mean_row_s1,
             output = out_s1)
strip_table_float(out_s1)
add_spanning_headers(out_s1)
cat("Saved:", out_s1, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S6: First-stage diagnostics -----
# Output: first_stage_diagnostics_neighbormun_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S6: First-stage diagnostics ===\n")

# Use the with-covariates models (indices 2, 4, 6, 8)
main_cov_models <- models_s1[c(2, 4, 6, 8)]
spec_names_s6 <- c("General population", "Falciparum API",
                    "Non-falciparum API", "Deaths per 1000")

out_s6 <- paste0(table_dir, "first_stage_diagnostics_neighbormun_area.tex")
make_iv_summary_table(main_cov_models, spec_names_s6, output_path = out_s6)
strip_table_float(out_s6)
cat("Saved:", out_s6, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S24: Indigenous population results -----
# Output: malaria_ip_neighbormun_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S24: Indigenous population ===\n")

ip_nocov <- feols(build_mu_iv_formula("malaria_indigenouspop22_api", "none"),
                  data = mod_df, cluster = "cd_mun")
ip_cov   <- feols(build_mu_iv_formula("malaria_indigenouspop22_api", "select"),
                  data = mod_df, cluster = "cd_mun")

mean_ip <- round(mean(mod_df$malaria_indigenouspop22_api, na.rm = TRUE), 2)
mean_row_ip <- data.frame(term = "Mean of DV", m1 = mean_ip, m2 = mean_ip)

out_s24 <- paste0(table_dir, "malaria_ip_neighbormun_area.tex")
msummary(list(ip_nocov, ip_cov),
         stars = TRUE,
         fmt = fmt_smart,
         coef_rename = coef_rename_mu,
         gof_map = gof_map_iv,
         add_rows = mean_row_ip,
         output = out_s24)
strip_table_float(out_s24)
cat("Saved:", out_s24, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S25: Garimpo results -----
# Output: malaria_garimpo_neighbormun_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S25: Garimpo ===\n")

garimpo_nocov <- feols(build_mu_iv_formula("n_garimpo", "none"),
                       data = mod_df, cluster = "cd_mun")
garimpo_cov   <- feols(build_mu_iv_formula("n_garimpo", "select"),
                       data = mod_df, cluster = "cd_mun")

mean_garimpo <- round(mean(mod_df$n_garimpo, na.rm = TRUE), 2)
mean_row_garimpo <- data.frame(term = "Mean of DV",
                               m1 = mean_garimpo, m2 = mean_garimpo)

out_s25 <- paste0(table_dir, "malaria_garimpo_neighbormun_area.tex")
msummary(list(garimpo_nocov, garimpo_cov),
         stars = TRUE,
         fmt = fmt_smart,
         coef_rename = coef_rename_mu,
         gof_map = gof_map_iv,
         add_rows = mean_row_garimpo,
         output = out_s25)
strip_table_float(out_s25)
cat("Saved:", out_s25, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S3: Health expenditure heterogeneity -----
# Output: malaria_general_hetero_effect_health_expenditure.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S3: Health expenditure heterogeneity ===\n")

hetero_outcomes <- c("malaria_allpop_api", "falciparum_api", "non_f_api",
                     "n_malaria_death_per1000")

# Panel A: without interaction (base models with health controls)
models_a <- lapply(hetero_outcomes, function(y) {
  fmla <- as.formula(paste0(
    y, " ~ ", mu_select_cov,
    " + health_expenditure_per_capita + hospital_visits_per1000",
    " | cd_mun + year | gold_mining_area ~ ", mu_instruments
  ))
  feols(fmla, data = mod_df, cluster = "cd_mun")
})
names(models_a) <- c("(1)", "(2)", "(3)", "(4)")

# Panel B: with health expenditure interaction
models_b <- lapply(hetero_outcomes, function(y) {
  fmla <- as.formula(paste0(
    y, " ~ ", mu_select_cov,
    " + health_expenditure_per_capita + hospital_visits_per1000",
    " | cd_mun + year",
    " | gold_mining_area + gold_mining_area:health_expenditure_per_capita",
    " ~ ", mu_instruments,
    " + pp3_delta_inc_2yrgp:health_expenditure_per_capita",
    " + pp3_gamma_mlp_2yrgp:health_expenditure_per_capita"
  ))
  feols(fmla, data = mod_df, cluster = "cd_mun")
})
names(models_b) <- c("(1)", "(2)", "(3)", "(4)")

# Coefficient display settings
coef_rename_health <- c(
  "fit_gold_mining_area" = "Gold mining area",
  "health_expenditure_per_capita" = "Health expenditure",
  "fit_gold_mining_area:health_expenditure_per_capita" = "Health expenditure * gold mining area",
  "hospital_visits_per1000" = "Outpatient visits per 1,000 pop."
)
coef_omit_health <- "forest_area_km|agriculture_area_km|pop|gdp_per_capita|mean_temp|annual_precip"

controls_row_a <- data.frame(
  term = "Environmental and socioeconomic controls",
  `(1)` = "X", `(2)` = "X", `(3)` = "X", `(4)` = "X",
  check.names = FALSE
)
attr(controls_row_a, "position") <- 7

controls_row_b <- data.frame(
  term = "Environmental and socioeconomic controls",
  `(1)` = "X", `(2)` = "X", `(3)` = "X", `(4)` = "X",
  check.names = FALSE
)
attr(controls_row_b, "position") <- 9

# Write each panel to temp file, extract body rows, combine
temp_a <- tempfile(fileext = ".tex")
temp_b <- tempfile(fileext = ".tex")

modelsummary(models_a, stars = TRUE, fmt = fmt_smart,
             coef_rename = coef_rename_health,
             coef_omit = coef_omit_health,
             add_rows = controls_row_a,
             gof_map = gof_map_iv, output = temp_a)
lines_a <- readLines(temp_a)
file.remove(temp_a)

modelsummary(models_b, stars = TRUE, fmt = fmt_smart,
             coef_rename = coef_rename_health,
             coef_omit = coef_omit_health,
             add_rows = controls_row_b,
             gof_map = gof_map_iv, output = temp_b)
lines_b <- readLines(temp_b)
file.remove(temp_b)

# Extract body rows (between column headers and \end{talltblr})
get_body_rows <- function(lines) {
  header_idx <- grep("^& \\(1\\)", lines)
  end_idx <- grep("\\\\end\\{talltblr\\}", lines)
  if (length(header_idx) == 0 || length(end_idx) == 0) return(lines)
  lines[(header_idx + 1):(end_idx - 1)]
}

body_a <- get_body_rows(lines_a)
body_b <- get_body_rows(lines_b)

gof_start_a <- grep("^Num\\.Obs", body_a)
gof_start_b <- grep("^Num\\.Obs", body_b)

coef_rows_a <- body_a[1:(gof_start_a - 1)]
gof_rows_a  <- body_a[gof_start_a:length(body_a)]
coef_rows_b <- body_b[1:(gof_start_b - 1)]
gof_rows_b  <- body_b[gof_start_b:length(body_b)]

n_coef_a <- length(coef_rows_a)
n_gof_a  <- length(gof_rows_a)
n_coef_b <- length(coef_rows_b)
n_gof_b  <- length(gof_rows_b)

combined_s3 <- c(
  "\\begin{table}",
  "\\centering",
  "\\begin{talltblr}[",
  "entry=none,label=none,",
  "note{}={+ p $< 0.1$, * p $< 0.05$, ** p $< 0.01$, *** p $< 0.001$},",
  "]",
  "{",
  "colspec={Q[]Q[]Q[]Q[]Q[]},",
  paste0("hline{1}={1-5}{solid, black, 0.1em},"),
  paste0("hline{2}={1-5}{solid, black, 0.05em},"),
  paste0("hline{", 1 + 1 + 1 + n_coef_a + n_gof_a + 1, "}={1-5}{solid, black, 0.05em},"),
  paste0("hline{", 1 + 1 + 1 + n_coef_a + n_gof_a + 1 + 1 + n_coef_b + n_gof_b + 1,
         "}={1-5}{solid, black, 0.1em},"),
  "column{2-5}={}{halign=c},",
  "column{1}={}{halign=l},",
  "}",
  "& (1) & (2) & (3) & (4) \\\\",
  paste0("\\SetCell[c=5]{l} \\textit{Panel A: Without interaction} & & & & \\\\"),
  "& Malaria API & Falciparum API & Non-falc. API & Death rate \\\\",
  coef_rows_a,
  gof_rows_a,
  paste0("\\SetCell[c=5]{l} \\textit{Panel B: With health expenditure interaction} & & & & \\\\"),
  "& Malaria API & Falciparum API & Non-falc. API & Death rate \\\\",
  coef_rows_b,
  gof_rows_b,
  "\\end{talltblr}",
  "\\end{table}"
)

out_s3 <- paste0(table_dir, "malaria_general_hetero_effect_health_expenditure.tex")
writeLines(combined_s3, out_s3)
strip_table_float(out_s3)
cat("Saved:", out_s3, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S8: Exclusion restriction — cross-sectional rock shares -----
# Output: exclusion_rockshare_2003_area_neighbor.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S8: Exclusion restriction (rock shares) ===\n")

df03 <- mod_df %>% filter(year == 2003)

excl_models_s8 <- list(
  "Percent agriculture" = lm(percent_ag ~ pp3_delta_inc + pp3_gamma_mlp, data = df03),
  "Percent soy"         = lm(percent_soy ~ pp3_delta_inc + pp3_gamma_mlp, data = df03),
  "Percent forest"      = lm(percent_forest ~ pp3_delta_inc + pp3_gamma_mlp, data = df03),
  "Avg. temp."          = lm(mean_temp ~ pp3_delta_inc + pp3_gamma_mlp, data = df03),
  "Precip."             = lm(annual_precip ~ pp3_delta_inc + pp3_gamma_mlp, data = df03)
)

out_s8 <- paste0(table_dir, "exclusion_rockshare_2003_area_neighbor.tex")
modelsummary(excl_models_s8, stars = TRUE,
             statistic = NULL,
             coef_rename = c(
               "pp3_delta_inc" = "Gabro Serra Comprida",
               "pp3_gamma_mlp" = "Granito Pepita"
             ),
             gof_map = gof_map_iv,
             output = out_s8)
strip_table_float(out_s8)
cat("Saved:", out_s8, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S9: Exclusion restriction — IV first-stage on covariates -----
# Output: exclusion_iv_firststage_area_neighbor.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S9: Exclusion restriction (IV first-stage) ===\n")

# IV models with covariates as endogenous variables
excl_covs <- c("agriculture_area_km", "soy_area_km", "forest_area_km",
               "mean_temp", "annual_precip")
excl_labels <- c("Agriculture", "Soy", "Forest", "Avg. temp.", "Precip.")

excl_iv_models <- lapply(excl_covs, function(cov) {
  fmla <- as.formula(paste0(
    "malaria_allpop_api ~ 1 | cd_mun + year | ", cov, " ~ ", mu_instruments
  ))
  feols(fmla, data = mod_df, cluster = "cd_mun")
})

# Extract first-stage summaries
excl_fs_summaries <- lapply(excl_iv_models, function(m) summary(m, stage = 1))
names(excl_fs_summaries) <- excl_labels

# Get first-stage F-stats from the IV models
f_vals_iv <- sapply(excl_iv_models, function(m) {
  sprintf("%.3f", fitstat(m, type = "ivf")[[1]]$stat)
})
names(f_vals_iv) <- excl_labels
f_rows_iv <- data.frame(c(list(" " = "First-stage F"), as.list(f_vals_iv)),
                        check.names = FALSE)

out_s9 <- paste0(table_dir, "exclusion_iv_firststage_area_neighbor.tex")
modelsummary(excl_fs_summaries, stars = TRUE,
             statistic = NULL,
             coef_rename = c(
               "pp3_delta_inc_2yrgp" = "Gabro Serra Comprida x 2 yr gold price",
               "pp3_gamma_mlp_2yrgp" = "Granito Pepita x 2 yr gold price"
             ),
             add_rows = f_rows_iv,
             gof_map = gof_map_iv,
             output = out_s9)
strip_table_float(out_s9)
cat("Saved:", out_s9, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S12: Rotemberg weights diagnostics -----
# Output: rotemberg_weights_diagnostics.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S12: Rotemberg weights ===\n")

# --- Compute Rotemberg weights from panel data ---
# Following the approach in iplc_r/3_analysis_mu/10_Rotemberge_weights.R,
# which uses bartik.weight::bw(). Here we use the pure-R implementation
# defined in 01_setup.R (compute_rotemberg_weights).

# Raw rock areas (time-invariant geology) and gold price per instrument-year
rock_vars <- c("pp3_delta_inc", "pp3_gamma_mlp")

# Reshape rock shares to wide: one column per (year x rock_type) pair
rockshare <- mod_df %>%
  select(year, cd_mun, all_of(rock_vars)) %>%
  distinct()
rockshare$year2 <- rockshare$year

rockshare_wide <- rockshare %>%
  pivot_wider(
    id_cols = c("cd_mun", "year"),
    names_from = "year2",
    values_from = all_of(rock_vars),
    names_glue = "t{year2}_{.value}",
    values_fill = 0
  ) %>%
  arrange(cd_mun, year)

# Gold price shifter: one value per (year x rock_type) pair
goldprice_exp <- tidyr::crossing(
  year = sort(unique(mod_df$year)),
  var  = rock_vars
) %>%
  left_join(mod_df %>% select(year, goldprice_2yr) %>% distinct(), by = "year") %>%
  arrange(year, var)

# Z column names in the order produced by pivot_wider (grouped by rock type).
# This matches the original bw() call which uses setdiff(names(rockshare_wide), index).
Z_cols <- setdiff(names(rockshare_wide), c("cd_mun", "year"))

stopifnot(length(Z_cols) == nrow(goldprice_exp))

# Sort mod_df to match rockshare_wide row order before residualizing
mod_df_sorted <- mod_df %>% arrange(cd_mun, year)

# Residualize y, x, and Z on municipality + year FE (FWL theorem)
y_res <- resid(feols(malaria_allpop_api ~ 0 | cd_mun + year, data = mod_df_sorted))
x_res <- resid(feols(gold_mining_area   ~ 0 | cd_mun + year, data = mod_df_sorted))

rockshare_wide_fac <- rockshare_wide %>%
  mutate(cd_mun = as.factor(cd_mun), year = as.factor(year))

Z_res <- sapply(Z_cols, function(zvar) {
  resid(feols(as.formula(paste0(zvar, " ~ 0 | cd_mun + year")),
              data = rockshare_wide_fac))
})

G_vec <- goldprice_exp$goldprice_2yr

# Compute weights
rw <- compute_rotemberg_weights(y = y_res, x = x_res, Z_mat = Z_res, G_vec = G_vec)

# Assemble into tibble with labels from goldprice_exp (matching bw() output)
bw <- tibble(
  year          = goldprice_exp$year,
  var           = goldprice_exp$var,
  goldprice_2yr = goldprice_exp$goldprice_2yr,
  alpha         = rw$alpha,
  beta          = rw$beta
)

# --- Verify against pre-computed CSV ---
bw_precomputed <- read_csv(paste0(inter_dir, "rotemberg_weights_bw_fe.csv"),
                           show_col_types = FALSE)
bw_check <- bw_precomputed %>%
  select(year, var, alpha_ref = alpha, beta_ref = beta) %>%
  left_join(bw %>% select(year, var, alpha, beta), by = c("year", "var"))
max_alpha_diff <- max(abs(bw_check$alpha - bw_check$alpha_ref))
max_beta_diff  <- max(abs(bw_check$beta  - bw_check$beta_ref))
cat(sprintf("  Verification vs pre-computed CSV: max |alpha diff| = %.2e, max |beta diff| = %.2e\n",
            max_alpha_diff, max_beta_diff))
stopifnot(max_alpha_diff < 1e-6, max_beta_diff < 1e-6)

# Summary statistics
n_pairs <- nrow(bw)
n_negative <- sum(bw$alpha < 0)
frac_negative <- n_negative / n_pairs
sum_positive <- sum(bw$alpha[bw$alpha > 0])
sum_negative <- sum(bw$alpha[bw$alpha < 0])

inst_shares <- bw %>%
  group_by(var) %>%
  summarise(alpha_sum = sum(alpha), .groups = "drop") %>%
  mutate(share = alpha_sum / sum(alpha_sum))

beta_overall <- sum(bw$alpha * bw$beta) / sum(bw$alpha)

# Instrument labels
instrument_labels <- c(
  "pp3_delta_inc" = "Gabro Serra Comprida $\\times$ gold price",
  "pp3_gamma_mlp" = "Granito Pepita $\\times$ gold price"
)

# Top 10 by |alpha|
top_k <- bw %>%
  mutate(abs_alpha = abs(alpha)) %>%
  arrange(desc(abs_alpha)) %>%
  head(10) %>%
  mutate(instrument = instrument_labels[var])

# Panel A rows
panel_a_rows <- c(
  sprintf("Number of year-instrument pairs & %d \\\\", n_pairs),
  sprintf("Negative weights (share of pairs) & %d / %d (%.0f\\%%) \\\\",
          n_negative, n_pairs, frac_negative * 100),
  sprintf("Sum of positive $\\hat{\\alpha}_k$ & %.3f \\\\", sum_positive),
  sprintf("Sum of negative $\\hat{\\alpha}_k$ & %.3f \\\\", sum_negative),
  sprintf("Overall IV estimate ($\\hat{\\beta}_{IV}$) & %.3f \\\\", beta_overall)
)
for (i in seq_len(nrow(inst_shares))) {
  inst_name <- instrument_labels[inst_shares$var[i]]
  panel_a_rows <- c(panel_a_rows,
    sprintf("Share from %s & %.3f \\\\", inst_name, inst_shares$share[i])
  )
}

# Panel B rows
panel_b_header <- "Year & Instrument & Gold price (\\$/oz) & $\\hat{\\alpha}_k$ & $\\hat{\\beta}_k$ \\\\"
panel_b_rows <- top_k %>%
  mutate(row = sprintf("%d & %s & %.0f & %.3f & %.3f \\\\",
                       year, instrument, goldprice_2yr, alpha, beta)) %>%
  pull(row)

# Assemble table
table_lines_s12 <- c(
  "\\begin{minipage}{\\textwidth}",
  "\\centering",
  "\\textbf{Panel A: Summary statistics}\\\\[0.5em]",
  "\\begin{tabular}{lc}",
  "\\toprule",
  panel_a_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{minipage}",
  "",
  "\\vspace{1em}",
  "",
  "\\begin{minipage}{\\textwidth}",
  "\\centering",
  "\\textbf{Panel B: Top 10 year-instrument pairs by $|\\hat{\\alpha}_k|$}\\\\[0.5em]",
  "\\begin{tabular}{clccc}",
  "\\toprule",
  panel_b_header,
  "\\midrule",
  panel_b_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{minipage}",
  "",
  "\\vspace{0.5em}",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  "\\textit{Notes:} Rotemberg weights ($\\hat{\\alpha}_k$) decompose the overall Bartik IV estimate into contributions from each year-instrument pair, following Goldsmith-Pinkham, Sorkin, and Swift (2020). The overall IV estimate equals $\\hat{\\beta}_{IV} = \\sum_k \\hat{\\alpha}_k \\hat{\\beta}_k$, where $\\hat{\\beta}_k$ is the just-identified IV estimate using only instrument $k$. Instruments are geological rock types (area in km\\textsuperscript{2}) interacted with the 2-year rolling average gold price. Municipality and year fixed effects are partialled out via the Frisch--Waugh--Lovell theorem.",
  "\\end{minipage}"
)

out_s12 <- paste0(table_dir, "rotemberg_weights_diagnostics.tex")
writeLines(table_lines_s12, out_s12)
strip_table_float(out_s12)
cat("Saved:", out_s12, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S13: Leave-one-out robustness -----
# Output: leave_one_out_robustness.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S13: Leave-one-out ===\n")

instruments_full   <- mu_instruments
instruments_dropIJ <- "pp3_gamma_m2_2yrgp + pp3_gamma_m3_2yrgp + pp3_gamma_m1_2yrgp"

drop_itaituba     <- 150360
drop_jacareacanga <- 150375

subsamples <- list(
  "Full sample"       = list(data = mod_df, iv = instruments_full),
  "Drop Itaituba"     = list(data = mod_df %>% filter(cd_mun != drop_itaituba),
                             iv = instruments_full),
  "Drop Jacareacanga" = list(data = mod_df %>% filter(cd_mun != drop_jacareacanga),
                             iv = instruments_full),
  "Drop both"         = list(data = mod_df %>% filter(!cd_mun %in% c(drop_itaituba, drop_jacareacanga)),
                             iv = instruments_dropIJ)
)

format_coef_loo <- function(coef, se, pval) {
  if (is.na(coef)) return("--")
  stars <- ifelse(pval < 0.001, "***",
           ifelse(pval < 0.01, "**",
           ifelse(pval < 0.05, "*",
           ifelse(pval < 0.1, "+", ""))))
  paste0(formatC(coef, format = "f", digits = 4), stars,
         " (", formatC(se, format = "f", digits = 4), ")")
}

loo_results <- list()
for (s_name in names(subsamples)) {
  df_sub <- subsamples[[s_name]]$data
  iv_str <- subsamples[[s_name]]$iv
  fmla <- as.formula(paste0(
    "malaria_allpop_api ~ 1 | cd_mun + year | gold_mining_area ~ ", iv_str
  ))
  m <- tryCatch(
    feols(fmla, data = df_sub, cluster = "cd_mun"),
    error = function(e) NULL
  )
  if (is.null(m)) {
    loo_results[[length(loo_results) + 1]] <- tibble(
      subsample = s_name, coefficient = NA_real_,
      std_error = NA_real_, p_value = NA_real_,
      first_stage_F = NA_real_, n_obs = nrow(df_sub)
    )
  } else {
    fs <- fitstat(m, type = "ivf")
    loo_results[[length(loo_results) + 1]] <- tibble(
      subsample = s_name,
      coefficient = coef(m)["fit_gold_mining_area"],
      std_error = se(m)["fit_gold_mining_area"],
      p_value = pvalue(m)["fit_gold_mining_area"],
      first_stage_F = fs[[1]]$stat,
      n_obs = nobs(m)
    )
  }
}
loo_df <- bind_rows(loo_results)

table_rows_loo <- loo_df %>%
  mutate(
    coef_str = mapply(format_coef_loo, coefficient, std_error, p_value),
    f_str = ifelse(is.na(first_stage_F), "--",
                   formatC(first_stage_F, format = "f", digits = 2))
  ) %>%
  select(subsample, coef_str, f_str, n_obs)
names(table_rows_loo) <- c("Sample", "Coef (SE)", "F-stat", "N")

kbl_loo <- kbl(table_rows_loo,
               format = "latex", booktabs = TRUE,
               label = "leave_one_out", escape = FALSE,
               align = c("l", rep("c", ncol(table_rows_loo) - 1))) %>%
  footnote(general = c(
    "IV estimates of gold mining area on health outcomes.",
    "Instruments are reselected using rLASSO when droping both Itaituba and Jacareacanga.",
    "No covariates. Municipality and year FE. Standard errors clustered at municipality level.",
    "$+$ p$<$0.1, * p$<$0.05, ** p$<$0.01, *** p$<$0.001."
  ), escape = FALSE, threeparttable = TRUE)

out_s13 <- paste0(table_dir, "leave_one_out_robustness.tex")
save_kable(kbl_loo, out_s13)
strip_table_float(out_s13)
cat("Saved:", out_s13, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S17: Population-weighted IV -----
# Output: malaria_general_weightedbypop_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S17: Population-weighted IV ===\n")

models_s17 <- list()
for (y in outcomes_main) {
  models_s17[[length(models_s17) + 1]] <- feols(
    build_mu_iv_formula(y, "none"),
    data = mod_df, weights = ~ pop, cluster = "cd_mun"
  )
  models_s17[[length(models_s17) + 1]] <- feols(
    build_mu_iv_formula(y, "select"),
    data = mod_df, weights = ~ pop, cluster = "cd_mun"
  )
}
names(models_s17) <- names(models_s1)

mean_row_s17 <- mean_row_s1  # Same mean values (unweighted means used)

out_s17 <- paste0(table_dir, "malaria_general_weightedbypop_area.tex")
modelsummary(models_s17, stars = TRUE,
             coef_rename = coef_rename_mu,
             gof_map = gof_map_iv,
             add_rows = mean_row_s17,
             output = out_s17)
strip_table_float(out_s17)
add_spanning_headers(out_s17)
cat("Saved:", out_s17, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S19: ManyIV robust estimators -----
# Output: manyiv_alloutcomes_neighbormun_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S19: ManyIV ===\n")

iv_formula_manyiv <- function(outcome) {
  as.formula(paste0(
    outcome, " ~ gold_mining_area + as.factor(cd_mun) + as.factor(year) | ",
    "pp3_delta_inc_2yrgp + pp3_gamma_mlp_2yrgp + as.factor(cd_mun) + as.factor(year)"
  ))
}

manyiv_outcomes <- c("malaria_allpop_api", "falciparum_api", "non_f_api",
                     "n_malaria_death_per1000")
manyiv_labels <- c("Malaria API", "Falciparum API",
                   "Non-falciparum API", "Malaria deaths per 1,000")

manyiv_results <- lapply(manyiv_outcomes, function(y) {
  message("Running ManyIV for: ", y)
  IVreg(iv_formula_manyiv(y), data = mod_df, inference = c("standard"))
})
names(manyiv_results) <- manyiv_outcomes

# State-clustered SE from fixest for TSLS
state_se_list <- sapply(manyiv_outcomes, function(y) {
  fmla <- as.formula(paste0(y, " ~ 1 | cd_mun + year | gold_mining_area ~ ", mu_instruments))
  fit <- feols(fmla, data = mod_df, cluster = "sigla_uf")
  se(fit)["fit_gold_mining_area"]
})
names(state_se_list) <- manyiv_outcomes

# Build combined table
estimator_types  <- c("tsls", "mbtsls", "liml")
estimator_labels <- c("TSLS", "MBTSLS", "LIML")

combined_manyiv <- list()
for (i in seq_along(manyiv_outcomes)) {
  est_df <- as.data.frame(manyiv_results[[manyiv_outcomes[i]]]$estimate)
  for (j in seq_along(estimator_types)) {
    if (estimator_types[j] %in% rownames(est_df)) {
      row <- est_df[estimator_types[j], , drop = FALSE]
      state_se_val <- if (estimator_types[j] == "tsls") state_se_list[[manyiv_outcomes[i]]] else NA
      combined_manyiv[[length(combined_manyiv) + 1]] <- data.frame(
        Estimator = estimator_labels[j],
        Estimate = row[, 1],
        `Conventional SE` = row[, 2],
        `Robust SE` = row[, 3],
        `HTE robust SE` = if (ncol(row) >= 4) row[, 4] else NA,
        `State-clustered SE` = state_se_val,
        row.names = NULL, check.names = FALSE
      )
    }
  }
}
tab_manyiv <- do.call(rbind, combined_manyiv)

fmt_manyiv <- function(x) {
  ifelse(is.na(x), "---",
         ifelse(abs(x) < 0.001 & x != 0, sprintf("%.2e", x),
                sprintf("%.3f", x)))
}
tab_manyiv[, 2:6] <- lapply(tab_manyiv[, 2:6], fmt_manyiv)

kbl_manyiv <- kableExtra::kbl(
  tab_manyiv, format = "latex", booktabs = TRUE,
  caption = "Many-IV Robust Estimates: General Population Malaria Outcomes. Conventional, robust, and HTE-robust SEs from ManyIV (not clustered). State-clustered SEs via fixest reported for TSLS only.",
  label = "manyiv",
  align = c("l", "r", "r", "r", "r", "r"),
  linesep = ""
) %>%
  kableExtra::kable_styling(latex_options = "hold_position") %>%
  kableExtra::pack_rows(index = setNames(
    rep(length(estimator_types), length(manyiv_labels)),
    manyiv_labels
  ))

out_s19 <- paste0(table_dir, "manyiv_alloutcomes_neighbormun_area.tex")
writeLines(as.character(kbl_manyiv), out_s19)
strip_table_float(out_s19)
cat("Saved:", out_s19, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S20: OLS estimates -----
# Output: ols_alloutcomes_neighbormun_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S20: OLS ===\n")

extract_ols_row <- function(model, spec_label) {
  coef_val <- coef(model)["gold_mining_area"]
  se_val <- se(model)["gold_mining_area"]
  data.frame(
    Specification = spec_label,
    Coefficient = coef_val,
    `Std. Error` = se_val,
    row.names = NULL, check.names = FALSE
  )
}

ols_outcome_models <- list(
  "Malaria API" = list(
    feols(malaria_allpop_api ~ gold_mining_area, data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(as.formula(paste0("malaria_allpop_api ~ gold_mining_area + ", mu_select_cov)),
          data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(malaria_allpop_api ~ gold_mining_area | cd_mun + year,
          data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(as.formula(paste0("malaria_allpop_api ~ gold_mining_area + ", mu_select_cov, " | cd_mun + year")),
          data = mod_df, vcov = vcov_cluster("cd_mun"))
  ),
  "Falciparum API" = list(
    feols(falciparum_api ~ gold_mining_area, data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(as.formula(paste0("falciparum_api ~ gold_mining_area + ", mu_select_cov)),
          data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(falciparum_api ~ gold_mining_area | cd_mun + year,
          data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(as.formula(paste0("falciparum_api ~ gold_mining_area + ", mu_select_cov, " | cd_mun + year")),
          data = mod_df, vcov = vcov_cluster("cd_mun"))
  ),
  "Non-falciparum API" = list(
    feols(non_f_api ~ gold_mining_area, data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(as.formula(paste0("non_f_api ~ gold_mining_area + ", mu_select_cov)),
          data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(non_f_api ~ gold_mining_area | cd_mun + year,
          data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(as.formula(paste0("non_f_api ~ gold_mining_area + ", mu_select_cov, " | cd_mun + year")),
          data = mod_df, vcov = vcov_cluster("cd_mun"))
  ),
  "Malaria death rate" = list(
    feols(n_malaria_death_per1000 ~ gold_mining_area, data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(as.formula(paste0("n_malaria_death_per1000 ~ gold_mining_area + ", mu_select_cov)),
          data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(n_malaria_death_per1000 ~ gold_mining_area | cd_mun + year,
          data = mod_df, vcov = vcov_cluster("cd_mun")),
    feols(as.formula(paste0("n_malaria_death_per1000 ~ gold_mining_area + ", mu_select_cov, " | cd_mun + year")),
          data = mod_df, vcov = vcov_cluster("cd_mun"))
  )
)

spec_labels_ols <- c("Naive OLS (no FE)", "Naive OLS + covariates",
                     "OLS + FE", "OLS + FE + covariates")

combined_ols <- list()
for (i in seq_along(ols_outcome_models)) {
  for (j in seq_along(ols_outcome_models[[i]])) {
    combined_ols[[length(combined_ols) + 1]] <- extract_ols_row(
      ols_outcome_models[[i]][[j]], spec_labels_ols[j]
    )
  }
}
tab_ols <- do.call(rbind, combined_ols)

fmt_coef_ols <- function(x) {
  x <- as.numeric(x)
  ifelse(abs(x) < 0.001 & x != 0, sprintf("%.6f", x), sprintf("%.3f", x))
}
fmt_se_ols <- function(x) {
  x <- as.numeric(x)
  ifelse(abs(x) < 0.001 & x != 0, sprintf("(%.6f)", x), sprintf("(%.3f)", x))
}
tab_ols$Coefficient <- fmt_coef_ols(tab_ols$Coefficient)
tab_ols$`Std. Error` <- fmt_se_ols(tab_ols$`Std. Error`)

pack_index_ols <- setNames(rep(4, length(names(ols_outcome_models))),
                           names(ols_outcome_models))

mean_ols <- c(
  mean(mod_df$malaria_allpop_api, na.rm = TRUE),
  mean(mod_df$falciparum_api, na.rm = TRUE),
  mean(mod_df$non_f_api, na.rm = TRUE),
  mean(mod_df$n_malaria_death_per1000, na.rm = TRUE)
)

kbl_ols <- kableExtra::kbl(
  tab_ols, format = "latex", booktabs = TRUE,
  caption = "OLS Estimates: General Population Malaria Outcomes",
  label = "ols_alloutcomes",
  align = c("l", "r", "r"),
  linesep = ""
) %>%
  kableExtra::kable_styling(latex_options = "hold_position") %>%
  kableExtra::pack_rows(index = pack_index_ols) %>%
  kableExtra::footnote(
    general = paste0(
      "Mean outcomes: Malaria API = ", sprintf("%.2f", mean_ols[1]),
      "; Falciparum API = ", sprintf("%.2f", mean_ols[2]),
      "; Non-falciparum API = ", sprintf("%.2f", mean_ols[3]),
      "; Malaria death rate = ", sprintf("%.4f", mean_ols[4]),
      ". Standard errors clustered at municipality level in parentheses."
    ),
    threeparttable = TRUE
  )

out_s20 <- paste0(table_dir, "ols_alloutcomes_neighbormun_area.tex")
writeLines(as.character(kbl_ols), out_s20)
strip_table_float(out_s20)
cat("Saved:", out_s20, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S22: Falsification tests (hospitalizations) -----
# Output: falsification_test_combined.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S22: Falsification (hospitalizations) ===\n")

hosp_df <- read_csv(paste0(data_dir, "hospitalization_by_municipality.csv"),
                    show_col_types = FALSE)

# Merge hospitalization data with main panel
mod_df_fals <- mod_df %>% left_join(hosp_df, by = c("cd_mun", "year"))

# Normalize to per 1,000 population
disease_cols <- c("Pneumonia", "Accidents", "Ulcer", "HIV",
                  "Tuberculosis", "Chronic_resp", "Dermatitis",
                  "Diabetes", "STD")

mod_df_fals <- mod_df_fals %>%
  mutate(across(all_of(intersect(disease_cols, colnames(mod_df_fals))),
                ~ replace_na(., 0))) %>%
  mutate(across(all_of(intersect(disease_cols, colnames(mod_df_fals))),
                ~ . / pop * 1000))

# Run IV for each disease (combined 9 diseases for table)
fals_diseases <- c("Chronic_resp", "Ulcer", "Tuberculosis", "Dermatitis",
                   "Accidents", "STD", "Pneumonia", "HIV", "Diabetes")
fals_labels <- c("Chr. Resp.", "Ulcer", "TB", "Derm.", "Accidents",
                 "STD", "Pneumonia", "HIV", "Diabetes")

fals_models <- lapply(fals_diseases, function(d) {
  fmla <- as.formula(paste0(d, " ~ 1 | cd_mun + year | gold_mining_area ~ ", mu_instruments))
  feols(fmla, data = mod_df_fals, cluster = "cd_mun")
})
names(fals_models) <- fals_labels

falsification_note <- c(
  "Outcomes are hospitalizations per 1,000 population.",
  "Chr. Resp. = Chronic Respiratory; TB = Tuberculosis; Derm. = Dermatitis; STD = sexually transmitted diseases."
)

out_s22 <- paste0(table_dir, "falsification_test_combined.tex")
msummary(fals_models, stars = TRUE,
         coef_rename = c(
           "fit_gold_mining_area" = "Gold mining area",
           "percent_lag_forest" = "Lagged percent forest in municipality",
           "mean_temp" = "Annual average temperature",
           "annual_precip" = "Annual precipitation"
         ),
         gof_map = gof_map_iv,
         notes = falsification_note,
         output = out_s22)
strip_table_float(out_s22)

# Fix siunitx negative sign rendering
tex_fals <- readLines(out_s22)
tex_fals <- gsub("\\\\num\\{-0\\.000\\}", "$-$0.000", tex_fals)
writeLines(tex_fals, out_s22)
cat("Saved:", out_s22, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S28: Health expenditure as IV outcome -----
# Output: iv_health_expenditure_outcome_neighbormun_area.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S28: Health expenditure as outcome ===\n")

iv_hepc <- feols(build_mu_iv_formula("health_expenditure_per_capita", "none"),
                 data = mod_df, cluster = "cd_mun")

mean_hepc <- round(mean(mod_df$health_expenditure_per_capita, na.rm = TRUE), 2)
mean_row_hepc <- data.frame(term = "Mean of DV", m1 = mean_hepc)

out_s28 <- paste0(table_dir, "iv_health_expenditure_outcome_neighbormun_area.tex")
modelsummary(list("Per capita health expenditure" = iv_hepc),
             stars = TRUE,
             coef_rename = coef_rename_mu,
             gof_map = gof_map_iv,
             add_rows = mean_row_hepc,
             output = out_s28)
strip_table_float(out_s28)
cat("Saved:", out_s28, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S4: Balance table (2003) -----
# Output: balance_table_2003.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S4: Balance table ===\n")

# Classify municipalities as gold-producing vs neighbor
gm <- mod_df %>%
  group_by(cd_mun) %>%
  summarise(total_gold = sum(gold_mining_area), .groups = "drop") %>%
  mutate(mun_status = ifelse(total_gold > 0, "gold_producing", "gold_neighbor")) %>%
  select(-total_gold)

df_bal <- mod_df %>% left_join(gm, by = "cd_mun")

# Balance variables
bal_vars <- c("pop", "pop_indigenous_2010", "area_km2",
              "percent_forest", "percent_soy", "percent_ag", "defor_rate",
              "health_expenditure", "health_expenditure_per_capita")

df_2003 <- df_bal %>% filter(year == 2003)
bal_vars_2003 <- bal_vars[sapply(bal_vars, function(v) !all(is.na(df_2003[[v]])))]

balancetbl <- df_2003 %>%
  mutate(mun_status = group_labels_mu[mun_status]) %>%
  select(mun_status, all_of(bal_vars_2003)) %>%
  tbl_summary(
    by = mun_status,
    label = as.list(var_labels_mu[intersect(names(var_labels_mu), bal_vars_2003)]),
    statistic = all_continuous() ~ "{mean} ({sd})",
    missing = "no"
  ) %>%
  add_p(test = all_continuous() ~ "t.test") %>%
  modify_header(
    label ~ "**Variable**",
    all_stat_cols() ~ "**{level}** (N = {n})"
  )

out_s4 <- paste0(table_dir, "balance_table_2003.tex")
balancetbl %>%
  as_gt() %>%
  gt::gtsave(filename = out_s4)
strip_table_float(out_s4)
cat("Saved:", out_s4, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Municipality analysis complete ===\n")
cat("Tables generated:\n")
for (f in list.files(table_dir, pattern = "\\.tex$")) {
  cat("  ", f, "\n")
}
