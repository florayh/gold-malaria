* =============================================================================
* Project: Amazon IPLC Mining & Malaria — Reproduction Package
* Purpose: Cluster LASSO instrument selection for clean control sample
*          (gold munis + non-gold munis >300km from all gold).
*          Writes selections to CSV.
*
* REQUIRES: Stata 17+, lassopack (ssc install lassopack)
*
* OUTPUT: data/intermediate/rlasso_selections_mu_clean_control.csv
* =============================================================================

* Run from reproduction/ directory: cd reproduction && stata-se -b do code/iv_selection/rlasso_clean_control.do
local data_dir "data/postfe"
local output_csv "data/intermediate/rlasso_selections_mu_clean_control.csv"

* =============================================================================
* IMPORT DATA
* =============================================================================
clear
capture confirm file "`data_dir'/mu_panel_clean_control_300km_postFE.csv"
if _rc {
    di "ERROR: PostFE file not found. Generate it first with R data prep."
    exit 601
}
import delimited "`data_dir'/mu_panel_clean_control_300km_postFE.csv"

* =============================================================================
* PANEL SETUP
* =============================================================================
xtset cd_mun year

tab year, gen(y_)
unab Y : y_*

* =============================================================================
* AUTO-DETECT INSTRUMENT VARIABLES
* =============================================================================
ds *_2yrgp
local potentialiv `r(varlist)'
local n_cand : word count `potentialiv'

di "=========================================="
di "Clean control sample (gold + >300km controls)"
di "Number of potential instruments: `n_cand'"
di "=========================================="

* =============================================================================
* RUN RLASSO
* =============================================================================
rlasso gold_mining_area `Y' `potentialiv', fe robust cluster(cd_mun) partial(`Y')

local sel = e(selected)
local n_sel : word count `sel'

di ""
di "SELECTED INSTRUMENTS (`n_sel'):"
di "`sel'"
di ""

* =============================================================================
* WRITE CSV OUTPUT
* =============================================================================
tempname fh
file open `fh' using "`output_csv'", write replace
file write `fh' "panel,n_candidates,n_selected,instruments" _n
file write `fh' "mu_clean_control,`n_cand',`n_sel',`sel'" _n
file close `fh'

di "Results saved to: `output_csv'"
