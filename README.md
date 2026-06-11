# Reproduction Package

**Paper:** "Gold mining increases malaria burden in the Brazilian Amazon"

---

## Overview

This package reproduces all tables, figures, and inline statistics from the main text and Supplementary Information.  

## Data Sources

The analysis panels were constructed from the sources described below. 

### Malaria Outcomes

Anonymized individual malaria case records from the **Malaria Epidemiological Surveillance System (SIVEP-Malária)**, accessed via Brazil's Access to Information request system (LAI). Cases were aggregated by municipality of infection, excluding imported cases and cure verification tests (LVC). The Annual Parasite Index (API) is calculated as cases per 1,000 population. Results are reported for all malaria, *P. falciparum*, and non-*P. falciparum* infections separately.

Malaria mortality data from the **Mortality Information System (SIM)**, identified by ICD-10 codes. Accessed via DATASUS.

### Indigenous Health Data

Indigenous malaria cases identified via the **SIVEP** locality classification (CATEGORI codes 12=Maloca, 29=Aldeia, 44=Área Indígena) and aggregated to Indigenous health subdistrict polo bases (DSEI polo bases). Annual population by polo base from the **Special Secretary of Indigenous Health (SESAI)**, accessed via Access to Information requests.

### Gold Mining Area

Annual gold mining area (km²) within municipalities and DSEI polo bases from **MapBiomas Brasil Collection 10** land use classification, based on Landsat Collection 2 Tier 1 imagery at 30m resolution. 

### Instruments

Geological rock formation areas from the **Brazilian Geological Service (CPRM-GeoSGB)**. The database contains 602 geological formation types; area (km²) of each formation was calculated by study unit. International gold price data from the **United States Geological Survey (USGS)** Historical Statistics for Mineral and Material Commodities. Instruments are constructed as rock_area × 2-year rolling average gold price for each formation type.

### Covariates

| Variable | Source |
|----------|--------|
| Population | IBGE yearly municipal estimates (via `brpop` R package) |
| GDP per capita | IBGE municipal GDP accounts |
| Precipitation | CHIRPS (Funk et al., 2015), processed in Google Earth Engine |
| Temperature | ERA5-Land (Muñoz-Sabater et al., 2021), processed in Google Earth Engine |
| Forest, agriculture, pasture area | MapBiomas Brasil Collection 10 |
| Deforestation rate | calculated based on MapBiomas |
| Health expenditure per capita | FINBRA municipal finance data (2004–2024) |
| Hospital visits per 1,000 pop | SIH-RD (Hospital Information System) |

### Falsification Outcomes

Hospitalization counts by disease classification, including chronic respiratory disease, ulcer, tuberculosis, dermatitis, accidents, STDs, pneumonia, HIV, and diabetes,  are accessed from the **Hospital Information System (SIH-RD)** via DATASUS.

DSEI health outcomes for falsification tests (child mortality by cause) from **SESAI** and **SINASC/SIM**, accessed via Access to Information requests and compiled at the polo base level.

### Distance to mining variables

**Indigenous village locations.** Point locations of Indigenous settlements (aldeias) from the **FUNAI GeoServer** (`aldeias_pontosPoint.shp`).

**General population.** Gridded population counts at 100m resolution from the **WorldPop** project (UN-adjusted, constrained national totals; file pattern `bra_pop_YYYY_CN_100m`). Population grids are available from 2015 onward; for years before 2015, the 2015 grid was used as a proxy.

**Distance calculation.** For each year:
- *Indigenous villages (Figure S1):* Euclidean distance from each aldeia point to the nearest gold mining pixel was computed using `terra::distance()`. An aldeia is classified as "within 5 km" if its distance to the nearest mining pixel is less than 5,000 meters.
- *General population (Figure S2):* The mining raster was aggregated from 30m to 100m resolution (to match the population grid), and distance-to-nearest-mining was computed for every 100m cell. Population within 5 km of mining was then calculated by summing the WorldPop pixel values for all cells within 5,000 meters of a mining pixel by municipality.

### Spatial Boundaries

Municipality boundaries from **IBGE**, simplified from the original shapefiles. DSEI polo base boundaries from **SESAI**. Both simplified using `st_simplify(dTolerance = 500)`.

## Software Requirements

### R (primary)

- **R 4.4+**
- Required packages: `tidyverse`, `fixest`, `modelsummary`, `kableExtra`, `gtsummary`, `janitor`, `ManyIV`, `patchwork`, `sf`, `ggnewscale`

```r
install.packages(c("tidyverse", "fixest", "modelsummary", "kableExtra",
                    "gtsummary", "janitor", "ManyIV", "patchwork",
                    "sf", "ggnewscale"))
```

### Stata (optional, for LASSO verification only)

- **Stata 17+** (SE or MP)
- Required package: `lassopack` (`ssc install lassopack`)

Stata scripts verify the cluster LASSO instrument selection (Belloni et al., 2016). The R analysis scripts hardcode the selected instruments, so Stata is not needed to reproduce results.

## Directory Structure

```
reproduction/
├── README.md               # This file
├── run_all.R                # Master script (sources scripts 02-07)
├── code/
│   ├── 01_setup.R           # Package loading, paths, helper functions
│   ├── 02_analysis_municipality.R   # Municipality-level analyses
│   ├── 03_analysis_dsei.R          # DSEI-level analyses
│   ├── 04_analysis_spillover.R     # Distance decay, clean control, spillover
│   ├── 05_analysis_extended_panel.R # 2003-2024 panel results
│   ├── 06_iv_sensitivity.R         # IV threshold sensitivity (Figure S5)
│   ├── 07_combined_outputs.R       # Table 1, all figures
│   └── iv_selection/               # IV selection: post-FE filtering + LASSO
│       ├── generate_postfe.R       # (Optional) Regenerate post-FE panels
│       ├── iv_sensitivity_postfe.R # (Optional) Figure S5 step 1: quantile sweep
│       ├── iv_sensitivity_rlasso.do# (Optional) Figure S5 step 2: rlasso sweep
│       ├── rlasso_municipality.do
│       ├── rlasso_dsei.do
│       ├── rlasso_municipality_extended.do
│       ├── rlasso_dsei_extended.do
│       ├── rlasso_distance_bands.do
│       ├── rlasso_clean_control.do
│       ├── rlasso_gold_direct_spillover.do
│       └── rlasso_dsei_falsification.do
├── data/
│   ├── panels/              # Analysis panels (CSV)
│   ├── spillover/           # SUTVA robustness panels (CSV)
│   ├── spatial/             # Simplified shapefiles for Figure 1 maps
│   ├── intermediate/        # Pre-computed intermediate data
│   └── postfe/              # Post-FE residualized panels (for Stata)
└── output/
    ├── tables/              # Generated LaTeX tables
    └── figures/             # Generated figures (PNG, PDF)
```

## How to Run

### Step 1: Run all analyses

```bash
cd reproduction
Rscript run_all.R
```

This sources scripts 02 through 07 sequentially. All outputs are written to `output/tables/` and `output/figures/`.

### Step 2 (optional): Verify LASSO instrument selection

From the `reproduction/` directory, run any Stata script:

```bash
cd reproduction
stata-se -b do code/iv_selection/rlasso_municipality.do
stata-se -b do code/iv_selection/rlasso_dsei.do
stata-se -b do code/iv_selection/rlasso_distance_bands.do
stata-se -b do code/iv_selection/rlasso_clean_control.do
stata-se -b do code/iv_selection/rlasso_gold_direct_spillover.do
stata-se -b do code/iv_selection/rlasso_dsei_falsification.do
```

These scripts are for verification only. Selected instruments are either hardcoded in the R scripts (`01_setup.R`, `05_analysis_extended_panel.R`) or available as pre-computed CSVs in `data/intermediate/`. 

### Step 3 (optional): Regenerate post-FE panels and IV sensitivity analysis

The `code/iv_selection/` directory contains scripts that show how the pre-computed data files were generated. Pre-computed versions are included, so these are optional and for verification only.

- `generate_postfe.R` (R): Regenerates the 18 post-FE residualized panels in `data/postfe/`.
- `iv_sensitivity_postfe.R` (R) + `iv_sensitivity_rlasso.do` (Stata): Regenerate `data/intermediate/rlasso_sweep_selections.csv` for Figure S5. 

## Output Map

### Main Text

| Output | File | Script |
|--------|------|--------|
| Table 1 | `combined_main_table.tex` | `07_combined_outputs.R` |
| Figure 1 | `figure1.png` | `07_combined_outputs.R` |
| Figure 2 | `figure2_coefficient_plot.png` | `07_combined_outputs.R` |
| Figure 3 | `distance_decay_malaria_api.pdf` | `04_analysis_spillover.R` |
| Figure 4 | `marginal_effect_health_expenditure.png` | `07_combined_outputs.R` |

### Supplementary Information

| Output | File | Script |
|--------|------|--------|
| Table S1 | `malaria_general_neighbormun_area.tex` | `02_analysis_municipality.R` |
| Table S2 | `malaria_dsei_neighbor_area.tex` | `03_analysis_dsei.R` |
| Table S3 | `malaria_general_hetero_effect_health_expenditure.tex` | `02_analysis_municipality.R` |
| Table S4 | `balance_table_2003.tex` | `02_analysis_municipality.R` |
| Table S5 | `balance_table_dsei_2003.tex` | `03_analysis_dsei.R` |
| Table S6 | `first_stage_diagnostics_neighbormun_area.tex` | `02_analysis_municipality.R` |
| Table S7 | `first_stage_diagnostics_dsei_neighbor_area.tex` | `03_analysis_dsei.R` |
| Table S8 | `exclusion_rockshare_2003_area_neighbor.tex` | `02_analysis_municipality.R` |
| Table S9 | `exclusion_iv_firststage_area_neighbor.tex` | `02_analysis_municipality.R` |
| Table S10 | `dsei_exclusion_rockshare_2003_area.tex` | `03_analysis_dsei.R` |
| Table S11 | `dsei_exclusion_iv_firststage.tex` | `03_analysis_dsei.R` |
| Table S12 | `rotemberg_weights_diagnostics.tex` | `02_analysis_municipality.R` |
| Table S13 | `leave_one_out_robustness.tex` | `02_analysis_municipality.R` |
| Table S14 | `table_distance_decay_cov.tex` | `04_analysis_spillover.R` |
| Table S15 | `table_gold_direct_spillover.tex` | `04_analysis_spillover.R` |
| Table S16 | `table_clean_control_general.tex` | `04_analysis_spillover.R` |
| Table S17 | `malaria_general_weightedbypop_area.tex` | `02_analysis_municipality.R` |
| Table S18 | `malaria_dsei_neighbor_weighted_area.tex` | `03_analysis_dsei.R` |
| Table S19 | `manyiv_alloutcomes_neighbormun_area.tex` | `02_analysis_municipality.R` |
| Table S20 | `ols_alloutcomes_neighbormun_area.tex` | `02_analysis_municipality.R` |
| Table S21 | `malaria_dsei_neighbor_area_drop2019.tex` | `03_analysis_dsei.R` |
| Table S22 | `falsification_test_combined.tex` | `02_analysis_municipality.R` |
| Table S23 | `dsei_neighbor_falsification.tex` | `03_analysis_dsei.R` |
| Table S24 | `malaria_ip_neighbormun_area.tex` | `02_analysis_municipality.R` |
| Table S25 | `malaria_garimpo_neighbormun_area.tex` | `02_analysis_municipality.R` |
| Table S26 | `malaria_03-24panel/malaria_api_compare_panels.tex` | `05_analysis_extended_panel.R` |
| Table S27 | `malaria_03-24panel/malaria_dsei_compare_panels.tex` | `05_analysis_extended_panel.R` |
| Table S28 | `iv_health_expenditure_outcome_neighbormun_area.tex` | `02_analysis_municipality.R` |
| Fig S1 | `aldeia_mining_exposure_trend_goldpolo.png` | `07_combined_outputs.R` |
| Fig S2 | `pop_mining_exposure_trend_goldmuni.png` | `07_combined_outputs.R` |
| Fig S3 | `gold_price.png` | `07_combined_outputs.R` |
| Fig S4 | `trends_goldarea_combined.png` | `07_combined_outputs.R` |
| Fig S5 | `fig_iv_threshold_sensitivity.pdf` | `06_iv_sensitivity.R` |
| Fig S6 | `sivep_race_completeness.png` | `07_combined_outputs.R` |

## Data Dictionary

### Main Analysis Panels (`data/panels/`)

| File | Description | Years |
|------|-------------|-------|
| `mu_panel_goldneighbor_0319.csv` | Municipality panel: gold-producing + neighbors | 2003-2019 |
| `dsei_panel_goldneighbor_0319.csv` | DSEI polo base panel: gold-producing + neighbors | 2003-2019 |
| `mu_panel_goldneighbor_full.csv` | Extended municipality panel | 2003-2024 |
| `dsei_panel_goldneighbor_full.csv` | Extended DSEI panel | 2003-2024 |
| `mu_panel_all_amazon.csv` | All Legal Amazon municipalities (trends only) | 2003-2019 |
| `dsei_panel_all_amazon.csv` | All DSEI polo bases (trends only) | 2003-2019 |
| `hospitalization_by_municipality.csv` | Hospital admission counts by disease (falsification) | 2003-2019 |
| `dsei_health_panel.csv` | DSEI health outcomes (falsification) | 2003-2019 |

### Spillover Panels (`data/spillover/`)

| File | Description |
|------|-------------|
| `mu_panel_band_{0_50,...,300_350}km.csv` | Distance band panels (7 bands) for Table S14 / Figure 3 |
| `mu_panel_clean_control_300km.csv` | Gold munis + controls >300km from gold municipalities (Table S16) |
| `mu_panel_gold_direct_spillover_200km.csv` | Direct + spillover joint estimation (Table S15) |

### Spatial Data (`data/spatial/`)

Simplified shapefiles (~3.4 MB total) for reproducing Figure 1 maps. Created from IBGE municipal boundaries and SESAI DSEI polo base boundaries using `st_simplify(dTolerance = 500)`.

| File | Description | Features |
|------|-------------|----------|
| `mun_amazon.gpkg` | Legal Amazon municipality boundaries | 808 municipalities |
| `polobase_amazon.gpkg` | DSEI polo base boundaries | 260 polo bases |
| `amazon_boundary.gpkg` | Legal Amazon state boundary outline | 1 polygon |

### Intermediate Data (`data/intermediate/`)

| File | Description |
|------|-------------|
| `rotemberg_weights_bw_fe.csv` | Pre-computed Rotemberg weights (Table S12) |
| `gold_price_usgs.csv` | USGS gold price series (Figure S3) |
| `aldeia_exposure_summary.csv` | Indigenous village mining exposure (Figure S1) |
| `pop_exposure_summary.csv` | Population mining exposure (Figure S2) |
| `rlasso_sweep_selections.csv` | Rlasso selections across quantile thresholds (Figure S5 input) |
| `iv_sensitivity_results.csv` | IV threshold sensitivity 2SLS results (Figure S5 output) |
| `sivep_race_completeness.png` | Pre-rendered SIVEP race completeness plot (Figure S6) |
| `rlasso_selections_*.csv` | LASSO-selected instruments for each specification |

### Post-FE Panels (`data/postfe/`)

Residualized panels for Stata LASSO verification. Each row is a unit-year observation after partialling out unit and year fixed effects. Used by scripts in `code/iv_selection/`. Can be regenerated from source panels using `code/iv_selection/generate_postfe.R`.

## Key Variables

### Treatment

| Variable | Description |
|----------|-------------|
| `gold_mining_area` | Gold mining area within municipality (km^2) |
| `goldmine_area` | Gold mining area near DSEI polo base (km^2) |
| `gold_mining_area_nb*` | Neighbor gold mining area at various distance bands |

### Outcomes

| Variable | Description |
|----------|-------------|
| `malaria_allpop_api` | Annual Parasite Index: malaria cases per 1,000 population |
| `falciparum_api` | Falciparum malaria API |
| `non_f_api` | Non-falciparum malaria API |
| `n_malaria_deaths` | Malaria deaths (raw count) |
| `malaria_api_by_polobase` | Malaria API at DSEI polo base level |




## Notes

- Figure 1 maps (panels A and B) are generated from simplified shapefiles included in `data/spatial/`. 
- Figure S6 is a pre-rendered PNG because it requires SIVEP malaria microdata (individual case records) not included in reproduction package.
