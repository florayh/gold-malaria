* =============================================================================
* Project: Amazon IPLC Mining & Malaria — Reproduction Package
* Purpose: Cluster LASSO instrument selection for direct + spillover model.
*          Two endogenous variables: own gold_mining_area and neighbor gold
*          (gold_mining_area_nb0_200km). Runs rlasso separately for each
*          using OWN vs NEIGHBOR IV pools. Writes selections to CSV.
*
* REQUIRES: Stata 17+, lassopack (ssc install lassopack)
*
* OUTPUT: data/intermediate/rlasso_selections_mu_spillover.csv
* =============================================================================

* Run from reproduction/ directory: cd reproduction && stata-se -b do code/iv_selection/rlasso_gold_direct_spillover.do
local data_dir "data/postfe"
local output_csv "data/intermediate/rlasso_selections_mu_spillover.csv"

* =============================================================================
* IMPORT DATA
* =============================================================================
clear
capture confirm file "`data_dir'/mu_panel_gold_direct_spillover_200km_postFE.csv"
if _rc {
    di "ERROR: PostFE file not found. Generate it first with R data prep."
    exit 601
}
import delimited "`data_dir'/mu_panel_gold_direct_spillover_200km_postFE.csv"

* =============================================================================
* PANEL SETUP
* =============================================================================
xtset cd_mun year

tab year, gen(y_)
unab Y : y_*

* =============================================================================
* AUTO-DETECT AND SEPARATE OWN vs NEIGHBOR IVs
* =============================================================================
ds *_2yrgp
local alliv `r(varlist)'

local owniv ""
local nbiv ""
foreach v of local alliv {
    if strpos("`v'", "_nb0_200km") > 0 {
        local nbiv `nbiv' `v'
    }
    else {
        local owniv `owniv' `v'
    }
}

local n_own : word count `owniv'
local n_nb : word count `nbiv'

di "=========================================="
di "Gold direct + spillover sample"
di "Own instruments: `n_own'"
di "Neighbor instruments: `n_nb'"
di "=========================================="

* =============================================================================
* Open CSV output
* =============================================================================
tempname fh
file open `fh' using "`output_csv'", write replace
file write `fh' "panel,n_candidates,n_selected,instruments" _n

* =============================================================================
* RLASSO 1: Own gold_mining_area (using OWN IVs)
* =============================================================================
di ""
di "RLASSO for OWN gold_mining_area (`n_own' instruments)"

rlasso gold_mining_area `Y' `owniv', fe robust cluster(cd_mun) partial(`Y')

local sel_own = e(selected)
local n_sel_own : word count `sel_own'

di "SELECTED (`n_sel_own'): `sel_own'"

file write `fh' "mu_spillover_own,`n_own',`n_sel_own',`sel_own'" _n

* =============================================================================
* RLASSO 2: Neighbor gold (using NEIGHBOR IVs)
* =============================================================================
di ""
di "RLASSO for NEIGHBOR gold_mining_area_nb0_200km (`n_nb' instruments)"

rlasso gold_mining_area_nb0_200km `Y' `nbiv', fe robust cluster(cd_mun) partial(`Y')

local sel_nb = e(selected)
local n_sel_nb : word count `sel_nb'

di "SELECTED (`n_sel_nb'): `sel_nb'"

file write `fh' "mu_spillover_nb,`n_nb',`n_sel_nb',`sel_nb'" _n

* =============================================================================
file close `fh'
di ""
di "============================================================"
di "Results saved to: `output_csv'"
di "============================================================"
