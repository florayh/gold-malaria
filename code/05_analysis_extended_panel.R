#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: Extended panel (2003-2024) analyses for municipality and DSEI.
#          Produces Tables S26, S27 (cross-panel comparison tables).
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/01_setup.R")

# Create subdirectory for extended panel tables
ext_table_dir <- paste0(table_dir, "malaria_03-24panel/")
dir.create(ext_table_dir, recursive = TRUE, showWarnings = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Municipality extended panel -----
# Instruments: pp3cc_2yrgp + mp3_delta_cs_2yrgp + pp3_delta_inc_2yrgp
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Municipality extended panel (03-24) ===\n")

mu_ext <- read_csv(paste0(data_dir, "mu_panel_goldneighbor_full.csv"),
                   show_col_types = FALSE)
mu_ext <- mu_ext %>%
  mutate(falciparum_api = 1000 * n_falciparum_mix / pop,
         non_f_api = 1000 * n_not_f / pop,
         n_malaria_death_per1000 = 1000 * n_malaria_deaths / pop)

cat("Loaded:", nrow(mu_ext), "rows,", n_distinct(mu_ext$cd_mun), "municipalities,",
    paste(range(mu_ext$year), collapse = "-"), "\n")

mu_ext_instruments <- "pp3cc_2yrgp + mp3_delta_cs_2yrgp + pp3_delta_inc_2yrgp"

# Version configs: full 03-24 and 03-24 excl COVID
mu_versions <- list(
  v0324 = list(name = "2003-2024", year_filter = 2003:2024, suffix = "_0324"),
  v0324_nc = list(name = "2003-2024 excl 2020-22",
                  year_filter = c(2003:2019, 2023:2024), suffix = "_0324nc")
)

mu_ext_results <- list()

for (vkey in names(mu_versions)) {
  v <- mu_versions[[vkey]]
  df_v <- mu_ext %>% filter(year %in% v$year_filter)

  cat("\n--- Municipality", v$name, ": N =", nrow(df_v), "---\n")

  make_fmla <- function(outcome, cov = NULL) {
    if (is.null(cov)) {
      as.formula(paste0(outcome, " ~ 1 | cd_mun + year | gold_mining_area ~ ", mu_ext_instruments))
    } else {
      as.formula(paste0(outcome, " ~ ", cov, " | cd_mun + year | gold_mining_area ~ ", mu_ext_instruments))
    }
  }

  # General population outcomes
  m_api_nocov  <- feols(make_fmla("malaria_allpop_api"), data = df_v, cluster = "cd_mun")
  m_api_cov    <- feols(make_fmla("malaria_allpop_api", mu_select_cov), data = df_v, cluster = "cd_mun")
  m_falc_nocov <- feols(make_fmla("falciparum_api"), data = df_v, cluster = "cd_mun")
  m_falc_cov   <- feols(make_fmla("falciparum_api", mu_select_cov), data = df_v, cluster = "cd_mun")
  m_nonf_nocov <- feols(make_fmla("non_f_api"), data = df_v, cluster = "cd_mun")
  m_nonf_cov   <- feols(make_fmla("non_f_api", mu_select_cov), data = df_v, cluster = "cd_mun")
  m_death_nocov <- feols(make_fmla("n_malaria_death_per1000"), data = df_v, cluster = "cd_mun")
  m_death_cov   <- feols(make_fmla("n_malaria_death_per1000", mu_select_cov), data = df_v, cluster = "cd_mun")

  gen_models <- list(
    "Malaria API" = m_api_nocov, "Malaria API" = m_api_cov,
    "Falciparum API" = m_falc_nocov, "Falciparum API" = m_falc_cov,
    "Non-falciparum API" = m_nonf_nocov, "Non-falciparum API" = m_nonf_cov,
    "Deaths per 1000" = m_death_nocov, "Deaths per 1000" = m_death_cov
  )

  mean_api   <- round(mean(df_v$malaria_allpop_api, na.rm = TRUE), 2)
  mean_falc  <- round(mean(df_v$falciparum_api, na.rm = TRUE), 2)
  mean_nonf  <- round(mean(df_v$non_f_api, na.rm = TRUE), 2)
  mean_death <- round(mean(df_v$n_malaria_death_per1000, na.rm = TRUE), 4)

  mean_row <- data.frame(term = "Mean of DV",
    m1 = mean_api, m2 = mean_api, m3 = mean_falc, m4 = mean_falc,
    m5 = mean_nonf, m6 = mean_nonf, m7 = mean_death, m8 = mean_death)

  mu_ext_results[[vkey]] <- list(api_nocov = m_api_nocov, api_cov = m_api_cov,
                                  mean_api = mean_api, n = nrow(df_v),
                                  fstat_api = fitstat(m_api_nocov, type = "ivf")[[1]]$stat)
}

# Table S26: Cross-panel comparison
if (length(mu_ext_results) >= 2) {
  compare_models <- list()
  compare_means <- c()
  for (vkey in names(mu_ext_results)) {
    v <- mu_versions[[vkey]]
    r <- mu_ext_results[[vkey]]
    compare_models[[v$name]] <- r$api_nocov
    compare_models[[paste0(v$name, " ")]] <- r$api_cov
    compare_means <- c(compare_means, r$mean_api, r$mean_api)
  }
  mean_row_cmp <- data.frame(term = "Mean of DV", t(compare_means))
  names(mean_row_cmp) <- c("term", paste0("m", seq_along(compare_means)))

  out_s26 <- paste0(ext_table_dir, "malaria_api_compare_panels.tex")
  msummary(compare_models, stars = TRUE, coef_rename = coef_rename_mu,
           gof_map = gof_map_iv, add_rows = mean_row_cmp, output = out_s26)
  strip_table_float(out_s26)
  add_spanning_headers(out_s26)
  cat("Saved:", out_s26, "\n")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DSEI extended panel -----
# Instruments: pp3cc_2yrgp + c2i_2yrgp + pp3_alfa_ar_2yrgp + pp4_gamma_po_2yrgp
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== DSEI extended panel (03-24) ===\n")

dsei_ext <- read_csv(paste0(data_dir, "dsei_panel_goldneighbor_full.csv"),
                     show_col_types = FALSE)
dsei_ext$seqid <- as.character(dsei_ext$seqid)

cat("Loaded:", nrow(dsei_ext), "rows,", n_distinct(dsei_ext$seqid), "polo bases,",
    paste(range(dsei_ext$year), collapse = "-"), "\n")

dsei_ext_instruments <- "pp3cc_2yrgp + c2i_2yrgp + pp3_alfa_ar_2yrgp + pp4_gamma_po_2yrgp"

dsei_versions <- list(
  v0324 = list(name = "2003-2024", year_filter = 2003:2024, suffix = "_0324"),
  v0324_nc = list(name = "2003-2024 excl 2020-22",
                  year_filter = c(2003:2019, 2023:2024), suffix = "_0324nc")
)

dsei_ext_results <- list()

for (vkey in names(dsei_versions)) {
  v <- dsei_versions[[vkey]]
  df_v <- dsei_ext %>% filter(year %in% v$year_filter)

  cat("\n--- DSEI", v$name, ": N =", nrow(df_v), "---\n")

  make_fmla_d <- function(outcome, cov = NULL) {
    if (is.null(cov)) {
      as.formula(paste0(outcome, " ~ 1 | seqid + year | goldmine_area ~ ", dsei_ext_instruments))
    } else {
      as.formula(paste0(outcome, " ~ ", cov, " | seqid + year | goldmine_area ~ ", dsei_ext_instruments))
    }
  }

  m_api_nocov <- feols(make_fmla_d("malaria_api_by_polobase"), data = df_v, cluster = "seqid")
  m_api_cov   <- feols(make_fmla_d("malaria_api_by_polobase", dsei_select_cov), data = df_v, cluster = "seqid")

  mean_api <- round(mean(df_v$malaria_api_by_polobase, na.rm = TRUE), 2)
  mean_row_d <- data.frame(term = "Mean of DV", m1 = mean_api, m2 = mean_api)

  dsei_models <- list("Malaria API" = m_api_nocov, "Malaria API" = m_api_cov)

  dsei_ext_results[[vkey]] <- list(api_nocov = m_api_nocov, api_cov = m_api_cov,
                                    mean_api = mean_api, n = nrow(df_v),
                                    fstat_api = fitstat(m_api_nocov, type = "ivf")[[1]]$stat)
}

# Table S27: Cross-panel comparison (DSEI)
if (length(dsei_ext_results) >= 2) {
  compare_models_d <- list()
  compare_means_d <- c()
  for (vkey in names(dsei_ext_results)) {
    v <- dsei_versions[[vkey]]
    r <- dsei_ext_results[[vkey]]
    compare_models_d[[v$name]] <- r$api_nocov
    compare_models_d[[paste0(v$name, " ")]] <- r$api_cov
    compare_means_d <- c(compare_means_d, r$mean_api, r$mean_api)
  }
  mean_row_cmp_d <- data.frame(term = "Mean of DV", t(compare_means_d))
  names(mean_row_cmp_d) <- c("term", paste0("m", seq_along(compare_means_d)))

  out_s27 <- paste0(ext_table_dir, "malaria_dsei_compare_panels.tex")
  msummary(compare_models_d, stars = TRUE, coef_rename = coef_rename_dsei,
           gof_map = gof_map_iv_dsei, add_rows = mean_row_cmp_d, output = out_s27)
  strip_table_float(out_s27)
  add_spanning_headers(out_s27)
  cat("Saved:", out_s27, "\n")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Extended panel analysis complete ===\n")
cat("Tables generated:\n")
cat("  S26: malaria_03-24panel/malaria_api_compare_panels.tex\n")
cat("  S27: malaria_03-24panel/malaria_dsei_compare_panels.tex\n")
