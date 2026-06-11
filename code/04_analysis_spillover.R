#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: SUTVA robustness analyses — distance decay, clean control, direct
#          spillover. Produces Tables S14, S15, S16 and Figure 3.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/01_setup.R")

# screen_instruments is defined in 01_setup.R

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S14 + Figure 3: Distance decay estimation -----
# Output: table_distance_decay_cov.tex, distance_decay_malaria_api.pdf
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S14 + Figure 3: Distance decay ===\n")

bands <- list(c(0,50), c(50,100), c(100,150), c(150,200),
              c(200,250), c(250,300), c(300,350))

# Read rlasso-selected instruments per band
rlasso_bands <- read_csv(paste0(inter_dir, "rlasso_selections_mu_bands.csv"),
                         show_col_types = FALSE)
instruments_by_band <- setNames(
  lapply(rlasso_bands$instruments, function(s) strsplit(s, " ")[[1]]),
  sub("^mu_band_", "", rlasso_bands$panel)
)

models_cov <- list()
band_stats <- list()

for (band in bands) {
  lo <- band[1]; hi <- band[2]
  band_label <- paste0(lo, "_", hi, "km")
  treatment <- paste0("gold_mining_area_nb", band_label)

  cat("\n--- Band:", band_label, "---\n")

  panel_path <- paste0(spill_dir, "mu_panel_band_", band_label, ".csv")
  if (!file.exists(panel_path)) { cat("  Panel not found. Skipping.\n"); next }
  mod_df <- read_csv(panel_path, show_col_types = FALSE)

  ivs <- instruments_by_band[[band_label]]
  if (is.null(ivs)) { cat("  No instruments found. Skipping.\n"); next }

  cat("  Screening", length(ivs), "instruments...\n")
  ivs_screened <- screen_instruments(mod_df, treatment, ivs)
  if (length(ivs_screened) == 0) { cat("  No instruments survived. Skipping.\n"); next }

  iv_str <- paste(ivs_screened, collapse = " + ")

  fmla_cov <- as.formula(paste0(
    "malaria_allpop_api ~ ", mu_select_cov, " | cd_mun + year | ",
    treatment, " ~ ", iv_str
  ))
  m_cov <- tryCatch(feols(fmla_cov, data = mod_df, cluster = "cd_mun"),
                     error = function(e) { cat("  Error:", e$message, "\n"); NULL })

  if (!is.null(m_cov)) {
    models_cov[[band_label]] <- m_cov
    fs <- fitstat(m_cov, type = "ivf")
    cat("  Coef:", round(coef(m_cov)[1], 3),
        " SE:", round(se(m_cov)[1], 3),
        " F:", round(fs[[1]]$stat, 1), "\n")

    band_stats[[band_label]] <- data.frame(
      band = band_label, midpoint = (lo + hi) / 2,
      coef = coef(m_cov)[1], se = se(m_cov)[1],
      n_mun = length(unique(mod_df$cd_mun)),
      n_obs = nobs(m_cov),
      f_stat = fs[[1]]$stat
    )
  }
}

# Table S14
if (length(models_cov) > 0) {
  coef_rename_bands <- setNames(
    rep("Neighbor gold mining (km\u00B2)", length(models_cov)),
    paste0("fit_gold_mining_area_nb", names(models_cov))
  )
  coef_rename_bands_cov <- c(coef_rename_bands,
    "forest_area_km" = "Forest area", "agriculture_area_km" = "Agriculture area",
    "pop" = "Population", "gdp_per_capita" = "GDP per capita",
    "mean_temp" = "Mean temperature", "annual_precip" = "Annual precipitation")

  out_s14 <- paste0(table_dir, "table_distance_decay_cov.tex")
  modelsummary(models_cov, stars = TRUE,
               coef_rename = coef_rename_bands_cov,
               gof_map = gof_map_iv,
               output = out_s14)
  strip_table_float(out_s14)
  cat("Saved:", out_s14, "\n")
}

# Figure 3: Distance decay plot
if (length(band_stats) > 0) {
  stats_df <- bind_rows(band_stats) %>%
    mutate(ci_lo = coef - 1.96 * se, ci_hi = coef + 1.96 * se,
           sig = ifelse(ci_lo > 0 | ci_hi < 0, "Significant", "Not significant"))

  p <- ggplot(stats_df, aes(x = midpoint, y = coef)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.8, linewidth = 0.6) +
    scale_x_continuous(
      breaks = stats_df$midpoint,
      labels = paste0(sub("_", "--", sub("km$", "", stats_df$band)), " km")
    ) +
    labs(x = "Distance band from focal municipality",
         y = "IV coefficient on neighbor gold mining area") +
    annotate("text", x = stats_df$midpoint, y = -0.5,
             label = paste0("n=", stats_df$n_mun), size = 3, color = "gray40") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 10))

  out_fig3 <- paste0(figure_dir, "distance_decay_malaria_api.pdf")
  ggsave(out_fig3, p, width = 8, height = 5)
  cat("Saved:", out_fig3, "\n")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S16: Clean control estimation -----
# Output: table_clean_control_general.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S16: Clean control ===\n")

# Read rlasso instruments for clean control
rlasso_cc <- read_csv(paste0(inter_dir, "rlasso_selections_mu_clean_control.csv"),
                      show_col_types = FALSE)
cc_ivs <- strsplit(rlasso_cc$instruments[1], " ")[[1]]
cat("Clean control rlasso instruments:", paste(cc_ivs, collapse = ", "), "\n")

cc_df <- read_csv(paste0(spill_dir, "mu_panel_clean_control_300km.csv"),
                  show_col_types = FALSE)
cc_df$falciparum_api <- 1000 * cc_df$n_falciparum_mix / cc_df$pop
cc_df$non_f_api <- 1000 * cc_df$n_not_f / cc_df$pop
cc_df$n_malaria_death_per1000 <- 1000 * cc_df$n_malaria_deaths / cc_df$pop

cat("Clean control sample:", nrow(cc_df), "obs,",
    length(unique(cc_df$cd_mun)), "municipalities\n")

# Screen instruments
cc_ivs_screened <- screen_instruments(cc_df, "gold_mining_area", cc_ivs)
cc_iv_str <- paste(cc_ivs_screened, collapse = " + ")

estimate_cc <- function(outcome) {
  fmla <- as.formula(paste0(outcome, " ~ ", mu_select_cov,
                            " | cd_mun + year | gold_mining_area ~ ", cc_iv_str))
  tryCatch(feols(fmla, data = cc_df, cluster = "cd_mun"),
           error = function(e) { cat("  Error:", e$message, "\n"); NULL })
}

res_malaria <- estimate_cc("malaria_allpop_api")
res_falc    <- estimate_cc("falciparum_api")
res_nonf    <- estimate_cc("non_f_api")
res_death   <- estimate_cc("n_malaria_death_per1000")

models_cc <- list(
  "Malaria API" = res_malaria,
  "Falciparum API" = res_falc,
  "Non-falciparum API" = res_nonf,
  "N deaths from malaria" = res_death
)
models_cc <- models_cc[!sapply(models_cc, is.null)]

if (length(models_cc) > 0) {
  mean_malaria_cc <- round(mean(cc_df$malaria_allpop_api, na.rm = TRUE), 2)
  mean_falc_cc    <- round(mean(cc_df$falciparum_api, na.rm = TRUE), 2)
  mean_nonf_cc    <- round(mean(cc_df$non_f_api, na.rm = TRUE), 2)
  mean_death_cc   <- round(mean(cc_df$n_malaria_death_per1000, na.rm = TRUE), 4)

  n_models <- length(models_cc)
  mean_vals <- sapply(names(models_cc), function(nm) {
    if (grepl("Falciparum", nm)) return(mean_falc_cc)
    if (grepl("Non-falc", nm)) return(mean_nonf_cc)
    if (grepl("death", nm, ignore.case = TRUE)) return(mean_death_cc)
    return(mean_malaria_cc)
  })
  mean_row_cc <- data.frame(term = "Mean of DV", t(mean_vals), check.names = FALSE)
  names(mean_row_cc) <- c("term", paste0("m", seq_len(n_models)))

  out_s16 <- paste0(table_dir, "table_clean_control_general.tex")
  modelsummary(models_cc, stars = TRUE, fmt = fmt_smart,
               coef_rename = coef_rename_mu,
               gof_map = gof_map_iv,
               add_rows = mean_row_cc,
               output = out_s16)
  strip_table_float(out_s16)
  add_spanning_headers(out_s16)
  cat("Saved:", out_s16, "\n")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table S15: Direct + spillover joint estimation -----
# Output: table_gold_direct_spillover.tex
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table S15: Direct + spillover ===\n")

# Read rlasso instruments (own and neighbor)
rlasso_spill <- read_csv(paste0(inter_dir, "rlasso_selections_mu_spillover.csv"),
                         show_col_types = FALSE)
rlasso_own_ivs <- strsplit(rlasso_spill$instruments[rlasso_spill$panel == "mu_spillover_own"], " ")[[1]]
rlasso_nb_ivs  <- strsplit(rlasso_spill$instruments[rlasso_spill$panel == "mu_spillover_nb"], " ")[[1]]

ds_df <- read_csv(paste0(spill_dir, "mu_panel_gold_direct_spillover_200km.csv"),
                  show_col_types = FALSE)
ds_df$falciparum_api <- 1000 * ds_df$n_falciparum_mix / ds_df$pop
ds_df$non_f_api <- 1000 * ds_df$n_not_f / ds_df$pop
ds_df$n_malaria_death_per1000 <- 1000 * ds_df$n_malaria_deaths / ds_df$pop

cat("Direct+spillover sample:", nrow(ds_df), "obs,",
    length(unique(ds_df$cd_mun)), "municipalities\n")

# Screen own and neighbor IVs separately
cat("Screening own IVs...\n")
own_ivs_screened <- screen_instruments(ds_df, "gold_mining_area", rlasso_own_ivs)
cat("Screening neighbor IVs...\n")
nb_ivs_screened <- screen_instruments(ds_df, "gold_mining_area_nb0_200km", rlasso_nb_ivs)

all_ivs <- c(own_ivs_screened, nb_ivs_screened)
iv_str_ds <- paste(all_ivs, collapse = " + ")

ds_outcomes <- c("malaria_allpop_api", "falciparum_api", "non_f_api", "n_malaria_death_per1000")
ds_labels   <- c("Malaria API", "Falciparum API", "Non-falciparum API", "Deaths per 1000")

models_ds <- list()
for (i in seq_along(ds_outcomes)) {
  fmla <- as.formula(paste0(
    ds_outcomes[i], " ~ ", mu_select_cov, " | cd_mun + year | ",
    "gold_mining_area + gold_mining_area_nb0_200km ~ ", iv_str_ds
  ))
  m <- tryCatch(feols(fmla, data = ds_df, cluster = "cd_mun"),
                error = function(e) { cat("  Error:", e$message, "\n"); NULL })
  if (!is.null(m)) models_ds[[ds_labels[i]]] <- m
}

if (length(models_ds) > 0) {
  coef_rename_ds <- c(
    "fit_gold_mining_area"           = "Own gold mining (km\u00B2)",
    "fit_gold_mining_area_nb0_200km" = "Neighbor gold mining 0--200km (km\u00B2)",
    "forest_area_km" = "Forest area (km\u00B2)", "agriculture_area_km" = "Agriculture area (km\u00B2)",
    "pop" = "Population", "gdp_per_capita" = "GDP per capita",
    "mean_temp" = "Mean temperature", "annual_precip" = "Annual precipitation"
  )

  n_models <- length(models_ds)
  mean_vals <- sapply(names(models_ds), function(nm) {
    if (grepl("Non-falc", nm)) return(round(mean(ds_df$non_f_api, na.rm = TRUE), 2))
    if (grepl("Falciparum", nm)) return(round(mean(ds_df$falciparum_api, na.rm = TRUE), 2))
    if (grepl("Death", nm)) return(round(mean(ds_df$n_malaria_death_per1000, na.rm = TRUE), 4))
    return(round(mean(ds_df$malaria_allpop_api, na.rm = TRUE), 2))
  })
  mean_row_ds <- data.frame(term = "Mean of DV", t(mean_vals), check.names = FALSE)
  names(mean_row_ds) <- c("term", paste0("m", seq_len(n_models)))

  out_s15 <- paste0(table_dir, "table_gold_direct_spillover.tex")
  modelsummary(models_ds, stars = TRUE, fmt = fmt_smart,
               coef_rename = coef_rename_ds,
               gof_map = gof_map_iv,
               add_rows = mean_row_ds,
               output = out_s15)
  strip_table_float(out_s15)
  add_spanning_headers(out_s15)
  cat("Saved:", out_s15, "\n")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Spillover analysis complete ===\n")
cat("Tables generated:\n")
cat("  S14: table_distance_decay_cov.tex\n")
cat("  S15: table_gold_direct_spillover.tex\n")
cat("  S16: table_clean_control_general.tex\n")
cat("Figures generated:\n")
cat("  Fig 3: distance_decay_malaria_api.pdf\n")
