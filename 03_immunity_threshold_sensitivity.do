*==================================================================*
* Immunity thresholds ≥10 / ≥100 / ≥1000 prevalence tables (svy)
* Outputs: $OUT/immunity_prevalence_thresholds.xlsx
*          long format: outcome group level_lab prev_pct l_pct u_pct
*==================================================================*
version 18
do "src/stata/_config.do"
use "$DERIVED/vk_evolve.dta", clear
svyset [pweight=totwt], singleunit(centered)

* --- Clean flags & numeric anti-HBs (if needed) ---
capture drop hbsab_num
gen double hbsab_num = .
replace hbsab_num = real(regexs(1)) if regexm(lower(hbsab_result),"([0-9.]+)")
gen byte hbsag_neg   = lower(hbsag_interpretation)=="nonreactive"
gen byte anti_hbc_pos = lower(anti_hbcii_interpretation)=="reactive"
label define LNO 0 "No" 1 "Yes", replace
label values hbsag_neg LNO
label values anti_hbc_pos LNO

* Keep rows with measurable anti-HBs for titres thresholds
drop if missing(hbsab_num)

* --- Threshold indicators (restrict to HBsAg-negative as per immunity defs) ---
foreach thr in 10 100 1000 {
    gen byte imm_any_`thr'      = (hbsag_neg==1 & hbsab_num>`thr')
    gen byte imm_vaccine_`thr'  = (hbsag_neg==1 & anti_hbc_pos==0 & hbsab_num>`thr')
    gen byte imm_resolved_`thr' = (hbsag_neg==1 & anti_hbc_pos==1 & hbsab_num>`thr')
    label values imm_any_`thr'      LNO
    label values imm_vaccine_`thr'  LNO
    label values imm_resolved_`thr' LNO
}

* --- Labelling for grouping variables ---
label define HBVDOB 0 "Born <1995" 1 "1995–1999" 2 "2000–2005", replace
capture label values hbvdob HBVDOB
label define HIV 0 "HIV-negative" 1 "PLWH", replace
capture label values hivstat HIV
label define SEX 1 "Male" 2 "Female", replace
capture label values sex SEX

* --- Helper: post results (svy mean + Wald CI) into long table ---
tempfile T
capture postclose P
postfile P str50 outcome str12 group str30 level_lab ///
    double(prev_pct l_pct u_pct) using "`T'", replace

program drop _all
program define __svy_post, rclass
    syntax varname, GROUPVAR(name) OUTCOME(string)
    tempvar z
    gen byte `z' = (`varlist'==1) if !missing(`varlist')
    quietly svy: mean `z'
    scalar b  = 100*_b[`z']
    scalar se = 100*_se[`z']
    scalar df = e(df_r)
    scalar t  = invttail(df, 0.025)
    scalar l  = b - abs(t)*se
    scalar u  = b + abs(t)*se
    return scalar prev = b
    return scalar l    = l
    return scalar u    = u
end

* --- Define outcomes and groups ---
local outcomes imm_any_10 imm_any_100 imm_any_1000 ///
               imm_vaccine_10 imm_vaccine_100 imm_vaccine_1000 ///
               imm_resolved_10 imm_resolved_100 imm_resolved_1000

local groups Overall sex hbvdob hivstat

* --- Overall (no stratifier) ---
foreach ov of local outcomes {
    quietly __svy_post `ov', groupvar(hbsag_neg) outcome("`ov'")
    
    if "`ov'" == "imm_any_10" local olabel "Any immunity >10"
    if "`ov'" == "imm_any_100" local olabel "Any immunity >100"
    if "`ov'" == "imm_any_1000" local olabel "Any immunity >1000"
    if "`ov'" == "imm_vaccine_10" local olabel "Vaccine-derived >10"
    if "`ov'" == "imm_vaccine_100" local olabel "Vaccine-derived >100"
    if "`ov'" == "imm_vaccine_1000" local olabel "Vaccine-derived >1000"
    if "`ov'" == "imm_resolved_10" local olabel "Resolved (infection-derived) >10"
    if "`ov'" == "imm_resolved_100" local olabel "Resolved (infection-derived) >100"
    if "`ov'" == "imm_resolved_1000" local olabel "Resolved (infection-derived) >1000"
    
    post P ("`olabel'") ("Overall") ("All") (r(prev)) (r(l)) (r(u))
}

* --- By sex, hbvdob, hivstat ---
foreach g of local groups {
    if "`g'"=="Overall" continue
    levelsof `g', local(levels)
    
    foreach ov of local outcomes {
        if "`ov'" == "imm_any_10" local olabel "Any immunity >10"
        if "`ov'" == "imm_any_100" local olabel "Any immunity >100"
        if "`ov'" == "imm_any_1000" local olabel "Any immunity >1000"
        if "`ov'" == "imm_vaccine_10" local olabel "Vaccine-derived >10"
        if "`ov'" == "imm_vaccine_100" local olabel "Vaccine-derived >100"
        if "`ov'" == "imm_vaccine_1000" local olabel "Vaccine-derived >1000"
        if "`ov'" == "imm_resolved_10" local olabel "Resolved (infection-derived) >10"
        if "`ov'" == "imm_resolved_100" local olabel "Resolved (infection-derived) >100"
        if "`ov'" == "imm_resolved_1000" local olabel "Resolved (infection-derived) >1000"
        
        foreach L of local levels {
            preserve
            keep if `g'==`L'
            quietly __svy_post `ov', groupvar(`g') outcome("`ov'")
            local lab : label (`g') `L'
            if "`lab'"=="" local lab "`L'"
            post P ("`olabel'") ("`g'") ("`lab'") (r(prev)) (r(l)) (r(u))
            restore
        }
    }
}

postclose P
use "`T'", clear
format prev_pct l_pct u_pct %6.1f
order outcome group level_lab prev_pct l_pct u_pct
sort outcome group level_lab

* --- Write to Excel (new file or replace) ---
export excel using "$OUT/immunity_prevalence_thresholds.xlsx", ///
    firstrow(variables) replace
di as res "Wrote: $OUT/immunity_prevalence_thresholds.xlsx"
