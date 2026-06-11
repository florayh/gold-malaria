* =============================================================================
* Project: Amazon IPLC Mining & Malaria — Reproduction Package
* Purpose: Cluster LASSO instrument selection for each distance band
*          (0-50, 50-100, ..., 300-350 km). Writes selections to CSV.
*
* REQUIRES: Stata 17+, lassopack (ssc install lassopack)
*
* OUTPUT: data/intermediate/rlasso_selections_mu_bands.csv
* =============================================================================

* Run from reproduction/ directory: cd reproduction && stata-se -b do code/iv_selection/rlasso_distance_bands.do
local data_dir "data/postfe"
local output_csv "data/intermediate/rlasso_selections_mu_bands.csv"

* Open output CSV
tempname fh
file open `fh' using "`output_csv'", write replace
file write `fh' "panel,n_candidates,n_selected,instruments" _n

* ==============================================================================
* Helper program
* ==============================================================================
capture program drop write_rlasso_result
program define write_rlasso_result
    args fh panel_name depvar potentialiv clustervar

    local n_cand : word count `potentialiv'

    tab year, gen(y_)
    unab Y : y_*

    rlasso `depvar' `Y' `potentialiv', fe robust cluster(`clustervar') partial(`Y')

    local sel = e(selected)
    local n_sel : word count `sel'

    di ""
    di "PANEL: `panel_name'"
    di "  Candidates: `n_cand'"
    di "  Selected (`n_sel'): `sel'"
    di ""

    file write `fh' "`panel_name',`n_cand',`n_sel',`sel'" _n

    drop y_*
end

* ==============================================================================
* Run all 7 bands
* ==============================================================================
foreach band in "0_50" "50_100" "100_150" "150_200" "200_250" "250_300" "300_350" {
    di "============================================================"
    di "Processing band: `band'km"
    di "============================================================"
    clear
    capture confirm file "`data_dir'/mu_panel_band_`band'km_postFE.csv"
    if _rc {
        di "  PostFE file not found for band `band'km — skipping"
        continue
    }
    import delimited "`data_dir'/mu_panel_band_`band'km_postFE.csv"
    xtset cd_mun year
    ds *_2yrgp
    local potentialiv `r(varlist)'
    write_rlasso_result `fh' "mu_band_`band'km" "gold_mining_area_nb`band'km" "`potentialiv'" "cd_mun"
}

* ==============================================================================
file close `fh'
di ""
di "============================================================"
di "ALL BANDS COMPLETE"
di "Results saved to: `output_csv'"
di "============================================================"
