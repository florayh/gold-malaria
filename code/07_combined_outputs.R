#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: Combined Table 1, all main-text and SI figures.
#          Produces Table 1, Figures 1-4, Figures S1-S6.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/01_setup.R")

# Additional packages for figures
library(patchwork)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Table 1: Combined main table (Panels A/B/C) -----
# Custom LaTeX — merges municipality IV, DSEI IV, and health expenditure
# heterogeneity results (with-covariates spec only)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Table 1: Combined main table ===\n")

# --- Load data ---
mu_df <- read_csv(paste0(data_dir, "mu_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
mu_df$falciparum_api <- 1000 * mu_df$n_falciparum_mix / mu_df$pop
mu_df$non_f_api <- 1000 * mu_df$n_not_f / mu_df$pop
mu_df$n_malaria_death_per1000 <- 1000 * mu_df$n_malaria_deaths / mu_df$pop

dsei_df <- read_csv(paste0(data_dir, "dsei_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
dsei_df$seqid <- as.character(dsei_df$seqid)
dsei_df$falciparum_api <- 1000 * dsei_df$malaria_n_falciparum_mix / dsei_df$population_polobase
dsei_df$non_f_api <- 1000 * dsei_df$malaria_n_not_f / dsei_df$population_polobase

# --- Helper functions ---
extract_coef <- function(model, coef_name) {
  ct <- coeftable(model)
  if (!coef_name %in% rownames(ct)) return(list(est = NA, se = NA, pval = NA))
  list(est = ct[coef_name, "Estimate"], se = ct[coef_name, "Std. Error"],
       pval = ct[coef_name, "Pr(>|t|)"])
}

format_stars <- function(pval) {
  if (is.na(pval)) return("")
  if (pval < 0.001) return("^{***}")
  if (pval < 0.01)  return("^{**}")
  if (pval < 0.05)  return("^{*}")
  if (pval < 0.1)   return("^{+}")
  return("")
}

fmt <- function(x, digits = 3) {
  if (is.na(x)) return("---")
  if (abs(x) > 0 & abs(x) < 0.005 & digits <= 3) digits <- 4
  if (abs(x) > 0 & abs(x) < 0.0005 & digits <= 4) digits <- 5
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

fmt_int <- function(x) {
  if (is.na(x)) return("---")
  formatC(x, format = "d", big.mark = ",")
}

get_fstat <- function(model) {
  fs <- fitstat(model, type = "ivf")
  round(fs[[1]]$stat, 1)
}

# --- Panel A: General population ---
mu_outcomes <- c("malaria_allpop_api", "falciparum_api", "non_f_api", "n_malaria_death_per1000")
panel_a <- lapply(mu_outcomes, function(y) {
  feols(build_mu_iv_formula(y, "select"), data = mu_df, cluster = "cd_mun")
})
names(panel_a) <- mu_outcomes

# --- Panel B: Indigenous population ---
dsei_outcomes <- c("malaria_api_by_polobase", "falciparum_api", "non_f_api")
panel_b <- lapply(dsei_outcomes, function(y) {
  feols(build_dsei_iv_formula(y, "select"), data = dsei_df, cluster = "seqid")
})
names(panel_b) <- dsei_outcomes

# --- Panel C: Health expenditure heterogeneity ---
panel_c <- lapply(mu_outcomes, function(y) {
  fmla <- as.formula(paste0(
    y, " ~ ", mu_select_cov,
    " + health_expenditure_per_capita + hospital_visits_per1000",
    " | cd_mun + year",
    " | gold_mining_area + gold_mining_area:health_expenditure_per_capita",
    " ~ ", mu_instruments,
    " + ", gsub("\\+", "+ ", gsub("(pp3_\\w+)", "\\1:health_expenditure_per_capita",
                                   mu_instruments))
  ))
  feols(fmla, data = mu_df, cluster = "cd_mun")
})
names(panel_c) <- mu_outcomes

# --- Extract coefficients ---
a_coefs <- lapply(panel_a, function(m) extract_coef(m, "fit_gold_mining_area"))
a_means <- sapply(mu_outcomes, function(y) mean(mu_df[[y]], na.rm = TRUE))

b_coefs <- lapply(panel_b, function(m) extract_coef(m, "fit_goldmine_area"))
b_means <- sapply(dsei_outcomes, function(y) mean(dsei_df[[y]], na.rm = TRUE))

c_coefs_main <- lapply(panel_c, function(m) extract_coef(m, "fit_gold_mining_area"))
c_coefs_int  <- lapply(panel_c, function(m) extract_coef(m, "fit_gold_mining_area:health_expenditure_per_capita"))

# --- Assemble LaTeX ---
lines <- c()
add <- function(...) lines <<- c(lines, paste0(...))

add("\\begin{table}[htbp]")
add("\\centering")
add("\\caption{Gold mining's effect on malaria incidence and mortality (2SLS-IV estimates)}")
add("\\label{tab:main_results}")
add("\\small")
add("\\begin{tabular}{l*{4}{c}}")
add("\\toprule")
add(" & (1) & (2) & (3) & (4) \\\\")
add(" & Malaria & \\textit{P.~falciparum} & Non-\\textit{falciparum} & Malaria deaths \\\\")
add(" & API & API & API & per 1,000 \\\\")
add("\\midrule")

# Panel A
add("\\multicolumn{5}{l}{\\textit{Panel A: General population (municipality-level)}} \\\\[3pt]")
add("Gold mining area (km\\textsuperscript{2}) & ",
    paste(sapply(1:4, function(i) {
      cc <- a_coefs[[i]]
      paste0("$", fmt(cc$est), format_stars(cc$pval), "$")
    }), collapse = " & "), " \\\\")
add(" & ",
    paste(sapply(1:4, function(i) paste0("(", fmt(a_coefs[[i]]$se), ")")), collapse = " & "), " \\\\[3pt]")
add("Observations & ",
    paste(sapply(panel_a, function(m) fmt_int(nobs(m))), collapse = " & "), " \\\\")
add("Mean of DV & ",
    paste(sapply(a_means, function(x) {
      if (abs(x) < 0.1) fmt(x, 3) else fmt(x, 2)
    }), collapse = " & "), " \\\\")
add("First-stage F & ",
    paste(sapply(panel_a, get_fstat), collapse = " & "), " \\\\")

add("\\midrule")

# Panel B
add("\\multicolumn{5}{l}{\\textit{Panel B: Indigenous population (indigenous health subdistrict-level)}} \\\\[3pt]")
b_est_row <- sapply(1:3, function(i) {
  cc <- b_coefs[[i]]
  paste0("$", fmt(cc$est), format_stars(cc$pval), "$")
})
add("Gold mining area (km\\textsuperscript{2}) & ",
    paste(c(b_est_row, "---"), collapse = " & "), " \\\\")
b_se_row <- sapply(1:3, function(i) paste0("(", fmt(b_coefs[[i]]$se), ")"))
add(" & ", paste(c(b_se_row, ""), collapse = " & "), " \\\\[3pt]")
add("Observations & ",
    paste(c(sapply(panel_b, function(m) fmt_int(nobs(m))), ""), collapse = " & "), " \\\\")
add("Mean of DV & ",
    paste(c(sapply(b_means, function(x) fmt(x, 2)), ""), collapse = " & "), " \\\\")
add("First-stage F & ",
    paste(c(sapply(panel_b, get_fstat), ""), collapse = " & "), " \\\\")

add("\\midrule")

# Panel C
add("\\multicolumn{5}{l}{\\textit{Panel C: Health expenditure heterogeneity (municipality-level)}} \\\\[3pt]")
add("Gold mining area (km\\textsuperscript{2}) & ",
    paste(sapply(1:4, function(i) {
      cc <- c_coefs_main[[i]]
      paste0("$", fmt(cc$est), format_stars(cc$pval), "$")
    }), collapse = " & "), " \\\\")
add(" & ",
    paste(sapply(1:4, function(i) paste0("(", fmt(c_coefs_main[[i]]$se), ")")), collapse = " & "), " \\\\[3pt]")
add("Health exp.~$\\times$ mining area & ",
    paste(sapply(1:4, function(i) {
      cc <- c_coefs_int[[i]]
      paste0("$", fmt(cc$est, 4), format_stars(cc$pval), "$")
    }), collapse = " & "), " \\\\")
add(" & ",
    paste(sapply(1:4, function(i) paste0("(", fmt(c_coefs_int[[i]]$se, 4), ")")), collapse = " & "), " \\\\[3pt]")
add("Observations & ",
    paste(sapply(panel_c, function(m) fmt_int(nobs(m))), collapse = " & "), " \\\\")
add("First-stage F & ",
    paste(sapply(panel_c, get_fstat), collapse = " & "), " \\\\")

add("\\midrule")
add("Unit \\& year FE & Yes & Yes & Yes & Yes \\\\")
add("Covariates & Yes & Yes & Yes & Yes \\\\")
add("Clustered SE & Yes & Yes & Yes & Yes \\\\")
add("\\bottomrule")
add("\\end{tabular}")
add("\\begin{tablenotes}[flushleft]")
add("\\footnotesize")
add("\\item \\textit{Notes:} Each cell reports a 2SLS-IV coefficient where gold mining area (km\\textsuperscript{2}) is instrumented by geological rock-type interactions with the rolling gold price. Panel~A estimates the effect on the general population at the municipality level (260 units, 2003--2019). Panel~B estimates the effect on indigenous populations at the indigenous health subdistrict level (106 units, 2003--2019); column~(4) is unavailable because indigenous-specific mortality data are not reported at this level. Panel~C adds an interaction between gold mining area and per-capita municipal health expenditure (R\\$) to test whether health system capacity mitigates mining's effect. All specifications include unit and year fixed effects, covariates (suppressed), and standard errors clustered at the unit level. $^{+}p<0.1$; $^{*}p<0.05$; $^{**}p<0.01$; $^{***}p<0.001$.")
add("\\end{tablenotes}")
add("\\end{table}")

out_table1 <- paste0(table_dir, "combined_main_table.tex")
writeLines(lines, out_table1)
cat("Saved:", out_table1, "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure 1: Maps + malaria trends composite -----
# Panel A: Municipality gold mining change map (2003-2019) from shapefiles
# Panel B: DSEI polo base gold mining change map (2003-2019) from shapefiles
# Panel C: Malaria API trends (gold vs non-gold, mun + DSEI)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure 1: Composite map + trends ===\n")

# Load full panels for trends and change computation
df_full_mun <- read_csv(paste0(data_dir, "mu_panel_all_amazon.csv"), show_col_types = FALSE)
df_full_mun$cd_mun <- as.character(df_full_mun$cd_mun)

df_full_dsei <- read_csv(paste0(data_dir, "dsei_panel_all_amazon.csv"), show_col_types = FALSE)

# Derive gold-producing vs neighbor IDs from the gold+neighbor panels
df_goldneighbor_mun <- read_csv(paste0(data_dir, "mu_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
df_goldneighbor_mun$cd_mun <- as.character(df_goldneighbor_mun$cd_mun)
gold_mun_ids <- df_goldneighbor_mun %>%
  group_by(cd_mun) %>% summarise(total = sum(gold_mining_area), .groups = "drop") %>%
  filter(total > 0) %>% pull(cd_mun)
neighbor_mun_ids <- setdiff(unique(df_goldneighbor_mun$cd_mun), gold_mun_ids)

df_goldneighbor_dsei <- read_csv(paste0(data_dir, "dsei_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
gold_polo_ids <- df_goldneighbor_dsei %>%
  group_by(seqid) %>% summarise(total = sum(goldmine_area), .groups = "drop") %>%
  filter(total > 0) %>% pull(seqid)
neighbor_polo_ids <- setdiff(unique(df_goldneighbor_dsei$seqid), gold_polo_ids)

# --- Load shapefiles ---
mun_amazon <- st_read(paste0(spatial_dir, "mun_amazon.gpkg"), quiet = TRUE)
polo_base_amz <- st_read(paste0(spatial_dir, "polobase_amazon.gpkg"), quiet = TRUE)
amazon_boundary <- st_read(paste0(spatial_dir, "amazon_boundary.gpkg"), quiet = TRUE)

# Classify municipalities: gold / neighbor / other
mun_amazon <- mun_amazon %>%
  mutate(mun_type = case_when(
    cd_mun_6d %in% gold_mun_ids ~ "gold",
    cd_mun_6d %in% neighbor_mun_ids ~ "neighbor",
    TRUE ~ "other"
  ))

other_mun <- mun_amazon %>% filter(mun_type == "other")
neighbor_mun <- mun_amazon %>% filter(mun_type == "neighbor")
gold_mun <- mun_amazon %>% filter(mun_type == "gold")

# Classify polo bases: gold / neighbor / other
polo_base_amz <- polo_base_amz %>%
  mutate(polo_type = case_when(
    seqid %in% gold_polo_ids ~ "gold",
    seqid %in% neighbor_polo_ids ~ "neighbor",
    TRUE ~ "other"
  ))

other_polo <- polo_base_amz %>% filter(polo_type == "other")
neighbor_polo <- polo_base_amz %>% filter(polo_type == "neighbor")
gold_polo <- polo_base_amz %>% filter(polo_type == "gold")

# --- Panel A: Municipality gold mining change map (2003-2019) ---
df_change_gold_mun <- df_full_mun %>%
  filter(cd_mun %in% gold_mun_ids, year %in% c(2003, 2019)) %>%
  select(cd_mun, year, gold_mining_area) %>%
  pivot_wider(names_from = year, values_from = gold_mining_area) %>%
  mutate(change_03_19 = `2019` - `2003`)

sf_change_gold_mun <- gold_mun %>%
  left_join(df_change_gold_mun, by = c("cd_mun_6d" = "cd_mun"))

max_abs_gold_mun <- max(abs(sf_change_gold_mun$change_03_19), na.rm = TRUE)

theme_change_map <- theme_void(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

p_map_mun <- ggplot() +
  geom_sf(data = other_mun, aes(fill = "Other Amazon"), color = "white", linewidth = 0.05) +
  geom_sf(data = neighbor_mun, aes(fill = "Neighbor"), color = "white", linewidth = 0.05) +
  scale_fill_manual(
    name = "", values = c("Other Amazon" = "grey90", "Neighbor" = "grey65"),
    guide = guide_legend(order = 2)
  ) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = sf_change_gold_mun, aes(fill = change_03_19), color = "white", linewidth = 0.1) +
  geom_sf(data = gold_mun, aes(color = "Gold-producing"), fill = NA, linewidth = 0.3,
          key_glyph = "path") +
  scale_fill_stepsn(
    colours = c("#92c5de", "#d1e5f0", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"),
    limits = c(-100, max_abs_gold_mun), n.breaks = 7,
    name = expression(paste("Change (", km^2, ")")),
    guide = guide_colorsteps(order = 3)
  ) +
  scale_color_manual(name = "", values = c("Gold-producing" = "grey30"),
                     guide = guide_legend(order = 1,
                                          override.aes = list(fill = NA, linewidth = 0.5, linetype = 1))) +
  geom_sf(data = amazon_boundary, fill = NA, color = "grey40", linewidth = 0.2) +
  coord_sf(datum = NA) +
  labs(title = "Gold mining area expansion by municipality: 2003-2019") +
  theme_change_map

# --- Panel B: DSEI polo base gold mining change map (2003-2019) ---
df_change_gold_dsei <- df_full_dsei %>%
  filter(seqid %in% gold_polo_ids, year %in% c(2003, 2019)) %>%
  select(seqid, year, goldmine_area) %>%
  pivot_wider(names_from = year, values_from = goldmine_area) %>%
  mutate(change_03_19 = `2019` - `2003`)

sf_change_gold_dsei <- gold_polo %>%
  left_join(df_change_gold_dsei, by = "seqid")

max_abs_gold_dsei <- max(abs(sf_change_gold_dsei$change_03_19), na.rm = TRUE)

p_map_dsei <- ggplot() +
  geom_sf(data = other_polo, aes(fill = "Other Amazon"), color = "white", linewidth = 0.05) +
  geom_sf(data = neighbor_polo, aes(fill = "Neighbor"), color = "white", linewidth = 0.05) +
  scale_fill_manual(
    name = "", values = c("Other Amazon" = "grey90", "Neighbor" = "grey65"),
    guide = guide_legend(order = 2)
  ) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = sf_change_gold_dsei, aes(fill = change_03_19), color = "white", linewidth = 0.1) +
  geom_sf(data = gold_polo, aes(color = "Gold-producing"), fill = NA, linewidth = 0.3,
          key_glyph = "path") +
  scale_fill_stepsn(
    colours = c("#92c5de", "#d1e5f0", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"),
    limits = c(-100, max_abs_gold_dsei), n.breaks = 7,
    name = expression(paste("Change (", km^2, ")")),
    guide = guide_colorsteps(order = 3)
  ) +
  scale_color_manual(name = "", values = c("Gold-producing" = "grey30"),
                     guide = guide_legend(order = 1,
                                          override.aes = list(fill = NA, linewidth = 0.5, linetype = 1))) +
  geom_sf(data = amazon_boundary, fill = NA, color = "grey40", linewidth = 0.2) +
  coord_sf(datum = NA) +
  labs(title = "Gold mining area expansion by indigenous health subdistrict: 2003-2019") +
  theme_change_map

# --- Panel C: Malaria trends ---
trends_mun <- df_full_mun %>%
  filter(year >= 2003, year <= 2019) %>%
  mutate(group = ifelse(cd_mun %in% gold_mun_ids,
                        "Gold-producing municipalities",
                        "Other municipalities")) %>%
  group_by(year, group) %>%
  summarise(mean_api = mean(malaria_allpop_api, na.rm = TRUE), .groups = "drop") %>%
  mutate(level = "Municipality")

trends_dsei <- df_full_dsei %>%
  filter(year >= 2003, year <= 2019) %>%
  mutate(group = ifelse(seqid %in% gold_polo_ids,
                        "Gold-producing indigenous health subdistricts",
                        "Other indigenous health subdistricts")) %>%
  group_by(year, group) %>%
  summarise(mean_api = mean(malaria_api_by_polobase, na.rm = TRUE), .groups = "drop") %>%
  mutate(level = "DSEI")

trends_all <- bind_rows(trends_mun, trends_dsei)

group_colors <- c(
  "Gold-producing municipalities" = "#d95f02",
  "Other municipalities" = "#d95f02",
  "Gold-producing indigenous health subdistricts" = "#1b9e77",
  "Other indigenous health subdistricts" = "#1b9e77"
)
group_linetypes <- c(
  "Gold-producing municipalities" = "solid",
  "Other municipalities" = "dashed",
  "Gold-producing indigenous health subdistricts" = "solid",
  "Other indigenous health subdistricts" = "dashed"
)
group_shapes <- c(
  "Gold-producing municipalities" = 16,
  "Other municipalities" = 1,
  "Gold-producing indigenous health subdistricts" = 17,
  "Other indigenous health subdistricts" = 2
)

p_trends <- ggplot(trends_all,
                   aes(x = year, y = mean_api, color = group,
                       linetype = group, shape = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(name = "", values = group_colors) +
  scale_linetype_manual(name = "", values = group_linetypes) +
  scale_shape_manual(name = "", values = group_shapes) +
  labs(x = "Year", y = "Malaria API (per 1,000)") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.margin = margin(t = -5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# --- Composite: maps on top, trends on bottom ---
fig1 <- (p_map_mun + p_map_dsei +
           plot_layout(widths = c(1, 1))) /
  p_trends +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C")),
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
  ) &
  theme(plot.tag.position = c(0, 1))

ggsave(paste0(figure_dir, "figure1.png"), fig1,
       width = 13, height = 12, dpi = 300)
cat("Saved:", paste0(figure_dir, "figure1.png"), "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure 2: Coefficient plot -----
# Parses LaTeX table outputs from Tables S1 and S2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure 2: Coefficient plot ===\n")

parse_talltblr_treatment <- function(tex_path) {
  lines <- readLines(tex_path)

  inner_close <- grep("tabularray inner close", lines)
  if (length(inner_close) == 0) stop("Cannot find 'tabularray inner close' in ", tex_path)

  header_idx <- NA
  for (i in (inner_close + 1):length(lines)) {
    if (grepl("&", lines[i])) { header_idx <- i; break }
  }
  if (is.na(header_idx)) stop("Cannot find header row in ", tex_path)

  coef_idx <- NA
  for (j in (header_idx + 1):length(lines)) {
    if (grepl("\\\\num\\{", lines[j])) { coef_idx <- j; break }
  }
  if (is.na(coef_idx)) stop("Cannot find coefficient row in ", tex_path)
  se_idx <- coef_idx + 1

  header_cells <- strsplit(lines[header_idx], "&")[[1]]
  header_cells <- trimws(gsub("\\\\\\\\", "", header_cells))
  header_cells <- trimws(gsub("\\\\SetCell\\[[^]]*\\]\\{[^}]*\\}", "", header_cells))
  outcome_names <- header_cells[header_cells != ""]

  coef_cells <- strsplit(lines[coef_idx], "&")[[1]]
  coef_cells <- trimws(gsub("\\\\\\\\", "", coef_cells))
  coef_values <- coef_cells[-1]

  if (length(outcome_names) * 2 == length(coef_values)) {
    outcomes <- rep(outcome_names, each = 2)
  } else {
    outcomes <- outcome_names
  }

  extract_coef_tex <- function(cell) {
    m <- regmatches(cell, regexec("\\\\num\\{([^}]+)\\}(\\*{0,3}\\+?)", cell))[[1]]
    if (length(m) < 2) return(list(estimate = NA_real_, stars = ""))
    list(estimate = as.numeric(m[2]), stars = if (length(m) >= 3) m[3] else "")
  }

  extract_se_tex <- function(cell) {
    m <- regmatches(cell, regexec("\\\\num\\{([^}]+)\\}", cell))[[1]]
    if (length(m) < 2) return(NA_real_)
    as.numeric(m[2])
  }

  se_cells <- strsplit(lines[se_idx], "&")[[1]]
  se_cells <- trimws(gsub("\\\\\\\\", "", se_cells))
  se_values <- se_cells[-1]

  n <- length(outcomes)
  tibble(
    outcome = outcomes,
    estimate = map_dbl(coef_values, ~ extract_coef_tex(.x)$estimate),
    stars = map_chr(coef_values, ~ extract_coef_tex(.x)$stars),
    se = map_dbl(se_values, extract_se_tex),
    conf.low = estimate - 1.96 * se,
    conf.high = estimate + 1.96 * se
  )
}

mun_coef_df <- parse_talltblr_treatment(paste0(table_dir, "malaria_general_neighbormun_area.tex")) %>%
  mutate(panel = "General population")

dsei_coef_df <- parse_talltblr_treatment(paste0(table_dir, "malaria_dsei_neighbor_area.tex")) %>%
  mutate(panel = "Indigenous health subdistrict")

# Label pairs (odd = no covariates, even = with covariates)
label_pairs <- function(df) {
  n <- nrow(df)
  labels <- character(n)
  for (i in seq(1, n, by = 2)) {
    base_name <- df$outcome[i]
    labels[i] <- base_name
    if (i + 1 <= n) labels[i + 1] <- paste0(base_name, " (cov)")
  }
  df$label <- labels
  df
}

mun_coef_df <- label_pairs(mun_coef_df)
dsei_coef_df <- label_pairs(dsei_coef_df)

mun_coef_df <- mun_coef_df %>%
  mutate(label = factor(label, levels = rev(c(
    "Malaria API", "Malaria API (cov)",
    "Non-falc. API", "Non-falc. API (cov)",
    "Falciparum API", "Falciparum API (cov)",
    "N deaths from malaria", "N deaths from malaria (cov)"
  ))))

dsei_coef_df <- dsei_coef_df %>%
  mutate(label = factor(label, levels = rev(c(
    "Malaria API", "Malaria API (cov)",
    "Non-falc. API", "Non-falc. API (cov)",
    "Falciparum API", "Falciparum API (cov)"
  ))))

color_map <- c(
  "Malaria API" = "#1b9e77", "Malaria API (cov)" = "#66c2a5",
  "Falciparum API" = "#d95f02", "Falciparum API (cov)" = "#fc8d62",
  "Non-falc. API" = "#7570b3", "Non-falc. API (cov)" = "#8da0cb",
  "N deaths from malaria" = "#e7298a", "N deaths from malaria (cov)" = "#f1a2c5"
)

coef_theme <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

all_data <- bind_rows(mun_coef_df, dsei_coef_df)
x_range <- c(min(all_data$conf.low, na.rm = TRUE), max(all_data$conf.high, na.rm = TRUE))
x_pad <- diff(x_range) * 0.05
x_limits <- c(x_range[1] - x_pad, x_range[2] + x_pad)

p_mun <- ggplot(mun_coef_df, aes(x = estimate, y = label, color = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0, linewidth = 0.6, orientation = "y") +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_map) +
  scale_x_continuous(limits = x_limits) +
  labs(title = "General population", x = NULL, y = NULL) +
  coef_theme

p_dsei <- ggplot(dsei_coef_df, aes(x = estimate, y = label, color = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0, linewidth = 0.6, orientation = "y") +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_map) +
  scale_x_continuous(limits = x_limits) +
  labs(title = "Indigenous health subdistrict", x = "Coefficient estimate and 95% CI", y = NULL) +
  coef_theme

p_fig2 <- p_mun / p_dsei +
  plot_layout(heights = c(4, 3)) +
  plot_annotation(
    title = "Effect of gold mining area on malaria outcomes",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(paste0(figure_dir, "figure2_coefficient_plot.png"),
       plot = p_fig2, width = 7, height = 7, dpi = 300)
cat("Saved:", paste0(figure_dir, "figure2_coefficient_plot.png"), "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure 3: Distance decay -----
# Already generated in 04_analysis_spillover.R — skip
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure 3: distance_decay_malaria_api.pdf ===\n")
if (file.exists(paste0(figure_dir, "distance_decay_malaria_api.pdf"))) {
  cat("  Already generated by 04_analysis_spillover.R\n")
} else {
  cat("  WARNING: Not found — run 04_analysis_spillover.R first\n")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure 4: Marginal effect of health expenditure -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure 4: Marginal effect of health expenditure ===\n")

# Re-run interaction model
fmla_health_cov3 <- as.formula(paste0(
  "malaria_allpop_api ~ ", mu_select_cov,
  " + health_expenditure_per_capita + hospital_visits_per1000|
  cd_mun + year|
  gold_mining_area  + gold_mining_area:health_expenditure_per_capita ~
  pp3_delta_inc_2yrgp + pp3_gamma_mlp_2yrgp +
  pp3_delta_inc_2yrgp:health_expenditure_per_capita + pp3_gamma_mlp_2yrgp:health_expenditure_per_capita"
))

health3 <- feols(fmla_health_cov3, data = mu_df, cluster = "cd_mun")

beta1 <- coef(health3)["fit_gold_mining_area"]
beta2 <- coef(health3)["fit_gold_mining_area:health_expenditure_per_capita"]

V <- vcov(health3)
v11 <- V["fit_gold_mining_area", "fit_gold_mining_area"]
v22 <- V["fit_gold_mining_area:health_expenditure_per_capita",
          "fit_gold_mining_area:health_expenditure_per_capita"]
v12 <- V["fit_gold_mining_area",
          "fit_gold_mining_area:health_expenditure_per_capita"]

W <- mu_df$health_expenditure_per_capita
w_grid <- seq(quantile(W, 0.01, na.rm = TRUE),
              quantile(W, 0.99, na.rm = TRUE),
              length.out = 200)

me <- beta1 + beta2 * w_grid
se_me <- sqrt(v11 + w_grid^2 * v22 + 2 * w_grid * v12)

me_df <- tibble(
  health_exp = w_grid,
  marginal_effect = me,
  ci_lower = me - 1.96 * se_me,
  ci_upper = me + 1.96 * se_me
)

p_top <- ggplot(me_df, aes(x = health_exp)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#d95f02", alpha = 0.2) +
  geom_line(aes(y = marginal_effect), color = "#d95f02", linewidth = 0.9) +
  labs(y = "Marginal effect of gold mining\non malaria API", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank())

hist_df <- tibble(health_exp = W) %>%
  filter(health_exp >= min(w_grid), health_exp <= max(w_grid))

p_bottom <- ggplot(hist_df, aes(x = health_exp)) +
  geom_histogram(bins = 40, fill = "grey70", color = "grey50", linewidth = 0.2) +
  labs(x = "Health expenditure per capita (R$)", y = "Frequency") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_fig4 <- p_top / p_bottom + plot_layout(heights = c(3, 1))

ggsave(paste0(figure_dir, "marginal_effect_health_expenditure.png"),
       p_fig4, width = 6, height = 5, dpi = 300)
cat("Saved:", paste0(figure_dir, "marginal_effect_health_expenditure.png"), "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure S1: Aldeia mining exposure trend -----
# Uses pre-computed summary CSV (GIS raster processing not reproducible here)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure S1: Aldeia mining exposure trend ===\n")

aldeia_summary_path <- paste0(inter_dir, "aldeia_exposure_summary.csv")
if (file.exists(aldeia_summary_path)) {
  aldeia_gold_summary <- read_csv(aldeia_summary_path, show_col_types = FALSE)

  p_s1 <- ggplot(aldeia_gold_summary, aes(x = year, y = pct_within_5km)) +
    geom_line(linewidth = 1.2, color = "#D55E00") +
    geom_point(size = 3, color = "#D55E00") +
    scale_x_continuous(breaks = seq(2003, 2024, by = 3)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = "Indigenous aldeias within 5km of gold mining",
      subtitle = "Gold-producing indigenous health subdistricts, 2003-2024",
      x = "Year",
      y = "% of aldeias within 5km of mining"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))

  ggsave(paste0(figure_dir, "aldeia_mining_exposure_trend_goldpolo.png"),
         p_s1, width = 8, height = 5, dpi = 300)
  cat("Saved:", paste0(figure_dir, "aldeia_mining_exposure_trend_goldpolo.png"), "\n")
} else {
  cat("  WARNING: aldeia_exposure_summary.csv not found — skipping Figure S1\n")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure S2: Population mining exposure trend -----
# Uses pre-computed summary CSV (GIS raster processing not reproducible here)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure S2: Population mining exposure trend ===\n")

pop_summary_path <- paste0(inter_dir, "pop_exposure_summary.csv")
if (file.exists(pop_summary_path)) {
  pop_gold_summary <- read_csv(pop_summary_path, show_col_types = FALSE)

  p_s2 <- ggplot(pop_gold_summary, aes(x = year, y = pct_within_5km)) +
    geom_line(linewidth = 1.2, color = "#0072B2") +
    geom_point(size = 3, color = "#0072B2") +
    scale_x_continuous(breaks = seq(2003, 2024, by = 3)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = "General population within 5km of gold mining",
      subtitle = "Gold-producing municipalities, 2003-2024",
      x = "Year",
      y = "% of population within 5km of mining"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))

  ggsave(paste0(figure_dir, "pop_mining_exposure_trend_goldmuni.png"),
         p_s2, width = 8, height = 5, dpi = 300)
  cat("Saved:", paste0(figure_dir, "pop_mining_exposure_trend_goldmuni.png"), "\n")
} else {
  cat("  WARNING: pop_exposure_summary.csv not found — skipping Figure S2\n")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure S3: Gold price trend -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure S3: Gold price trend ===\n")

goldprice <- read_csv(paste0(inter_dir, "gold_price_usgs.csv"), show_col_types = FALSE) %>%
  clean_names() %>%
  filter(year >= 2000)

p_s3 <- ggplot(goldprice, aes(x = year, y = gold_troy_ounce_price)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = "Year", y = "Gold Price (USD per troy ounce)",
       title = "International Gold Price Over Time") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(paste0(figure_dir, "gold_price.png"), p_s3, width = 8, height = 5, dpi = 300)
cat("Saved:", paste0(figure_dir, "gold_price.png"), "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure S4: Gold mining area trends (mun + DSEI combined) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure S4: Gold area trends ===\n")

goldarea_mun <- df_full_mun %>%
  group_by(year) %>%
  summarise(total_gold_area = sum(gold_mining_area, na.rm = TRUE), .groups = "drop") %>%
  mutate(level = "Municipalities")

goldarea_dsei <- df_full_dsei %>%
  group_by(year) %>%
  summarise(total_gold_area = sum(goldmine_area, na.rm = TRUE), .groups = "drop") %>%
  mutate(level = "Indigenous health subdistricts")

goldarea_all <- bind_rows(goldarea_mun, goldarea_dsei)

p_s4 <- ggplot(goldarea_all,
               aes(x = year, y = total_gold_area, color = level, shape = level)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(
    name = "",
    values = c("Municipalities" = "#d95f02",
               "Indigenous health subdistricts" = "#1b9e77")
  ) +
  scale_shape_manual(
    name = "",
    values = c("Municipalities" = 16,
               "Indigenous health subdistricts" = 17)
  ) +
  labs(x = "Year",
       y = expression(paste("Total gold mining area (", km^2, ")"))) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.margin = margin(t = -5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

ggsave(paste0(figure_dir, "trends_goldarea_combined.png"), p_s4,
       width = 8, height = 5, dpi = 300)
cat("Saved:", paste0(figure_dir, "trends_goldarea_combined.png"), "\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure S6: SIVEP race completeness -----
# Pre-rendered figure (requires raw SIVEP microdata not in reproduction package)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure S6: SIVEP race completeness ===\n")

sivep_src <- paste0(inter_dir, "sivep_race_completeness.png")
if (file.exists(sivep_src)) {
  file.copy(sivep_src, paste0(figure_dir, "sivep_race_completeness.png"), overwrite = TRUE)
  cat("Copied pre-rendered:", paste0(figure_dir, "sivep_race_completeness.png"), "\n")
} else {
  cat("  WARNING: sivep_race_completeness.png not found in intermediate data\n")
  cat("  This figure requires raw SIVEP microdata not included in the reproduction package.\n")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Combined outputs complete ===\n")
cat("Tables generated:\n")
cat("  Table 1: combined_main_table.tex\n")
cat("Figures generated:\n")
cat("  Fig 1:  figure1.png\n")
cat("  Fig 2:  figure2_coefficient_plot.png\n")
cat("  Fig 3:  distance_decay_malaria_api.pdf (from 04_analysis_spillover.R)\n")
cat("  Fig 4:  marginal_effect_health_expenditure.png\n")
cat("  Fig S1: aldeia_mining_exposure_trend_goldpolo.png\n")
cat("  Fig S2: pop_mining_exposure_trend_goldmuni.png\n")
cat("  Fig S3: gold_price.png\n")
cat("  Fig S4: trends_goldarea_combined.png\n")
cat("  Fig S5: fig_iv_threshold_sensitivity.pdf (from 06_iv_sensitivity.R)\n")
cat("  Fig S6: sivep_race_completeness.png\n")
