* =============================================================================
* Project: Amazon IPLC Mining & Malaria -- Reproduction Package
* Purpose: Cluster LASSO instrument selection for the DSEI falsification
*          mortality subsample (births >= 5, used for Table S23).
*          Writes selections to CSV for R to read directly.
*
* REQUIRES: Stata 17+, lassopack (ssc install lassopack)
*
* OUTPUT: data/intermediate/rlasso_selections_dsei_falsification.csv
* =============================================================================

* Run from reproduction/ directory: cd reproduction && stata-se -b do code/iv_selection/rlasso_dsei_falsification.do
local data_dir "data/postfe"
local output_csv "data/intermediate/rlasso_selections_dsei_falsification.csv"

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
    di "  Observations: `e(N)'"
    di "  Candidates: `n_cand'"
    di "  Selected (`n_sel'): `sel'"
    di ""

    file write `fh' "`panel_name',`n_cand',`n_sel',`sel'" _n

    drop y_*
end

* ==============================================================================
* Run falsification subsample
* ==============================================================================
foreach sample in "mortality" {
    di "============================================================"
    di "Processing subsample: `sample'"
    di "============================================================"
    clear
    capture confirm file "`data_dir'/dsei_falsification_`sample'_postFE.csv"
    if _rc {
        di "  PostFE file not found for `sample' -- skipping"
        continue
    }
    import delimited "`data_dir'/dsei_falsification_`sample'_postFE.csv"
    xtset seqid year
    ds *_2yrgp
    local potentialiv `r(varlist)'
    write_rlasso_result `fh' "`sample'" "goldmine_area" "`potentialiv'" "seqid"
}

* ==============================================================================
file close `fh'
di ""
di "============================================================"
di "FALSIFICATION SUBSAMPLE COMPLETE"
di "Results saved to: `output_csv'"
di "============================================================"
