/******************************************************************
 PI schematic: classify immunity/infection from HBsAg, anti-HBs, anti-HBc
******************************************************************/

version 18
do "src/stata/_config.do"
use "$DERIVED/vk_evolve.dta", clear
svyset [pweight=totwt], singleunit(centered)

tempvar hbsag_i hbsab_i ahbc_i
gen `hbsag_i' = lower(trim(hbsag_interpretation))
gen `hbsab_i' = lower(trim(hbsab_interpretation))
gen `ahbc_i'  = lower(trim(anti_hbcii_interpretation))

tempvar miss_hbsag miss_hbsab miss_ahbc
gen byte `miss_hbsag' = inlist(`hbsag_i',"","insufficient","equivocal")
gen byte `miss_hbsab' = inlist(`hbsab_i',"","insufficient","equivocal")
gen byte `miss_ahbc'  = inlist(`ahbc_i' ,"","insufficient","equivocal")

label define yesno 0 "No" 1 "Yes", replace

* numeric S/CO if available
capture confirm numeric variable hbv_num
if _rc {
    capture confirm string variable hbsag_result
    if !_rc gen double hbv_num = real(regexs(1)) if regexm(hbsag_result,"([0-9.]+)")
}

/******************************************************************
  B) Antibody flags (reused by both thresholds)
******************************************************************/
cap drop anti_hbs_pos anti_hbc_pos
gen byte anti_hbs_pos = (`hbsab_i'=="reactive") if !`miss_hbsab'
gen byte anti_hbc_pos = (`ahbc_i' =="reactive") if !`miss_ahbc'
label values anti_hbs_pos yesno
label values anti_hbc_pos yesno

/******************************************************************
  C) HBsAg positivity flags for BOTH thresholds
     - >1  (matches your current workflow)
     - >10 (new)
     Fallback to qualitative interpretation where S/CO missing
******************************************************************/
* >1
cap drop hbsag_gt1
gen byte hbsag_gt1 = .
replace hbsag_gt1 = (hbv_num > 1) if !missing(hbv_num)
replace hbsag_gt1 = 1 if missing(hbsag_gt1) & `hbsag_i'=="reactive"    & !`miss_hbsag'
replace hbsag_gt1 = 0 if missing(hbsag_gt1) & `hbsag_i'=="nonreactive" & !`miss_hbsag'
label values hbsag_gt1 yesno
label var     hbsag_gt1 "HBsAg positive (S/CO > 1 or interpretation)"

* >10
cap drop hbsag_gt10
gen byte hbsag_gt10 = .
replace hbsag_gt10 = (hbv_num > 10) if !missing(hbv_num)
replace hbsag_gt10 = 1 if missing(hbsag_gt10) & `hbsag_i'=="reactive"    & !`miss_hbsag'
replace hbsag_gt10 = 0 if missing(hbsag_gt10) & `hbsag_i'=="nonreactive" & !`miss_hbsag'
label values hbsag_gt10 yesno
label var     hbsag_gt10 "HBsAg positive (S/CO > 10 or interpretation)"

/******************************************************************
  D) PI schematic categories for BOTH thresholds
******************************************************************/
label define piimm 1 "Living with HBV infection (HBsAg+)" ///
                   2 "Immunised (never infected): HBsAg-/anti-HBs+/anti-HBc-" ///
                   3 "Cleared infection (± immunised): HBsAg-/anti-HBs+/anti-HBc+" ///
                   4 "Current or previous infection: HBsAg-/anti-HBs-/anti-HBc+" ///
                   5 "Never exposed & never immunised: HBsAg-/anti-HBs-/anti-HBc-", replace

* helper program to build the 5-way category from an HBsAg flag
capture program drop __make_pi_cat
program define __make_pi_cat
    syntax varname, OUT(name)
    cap drop `out'
    gen byte `out' = .
    replace `out' = 1 if `varlist'==1
    replace `out' = 2 if missing(`out') & `varlist'==0 & anti_hbs_pos==1 & anti_hbc_pos==0
    replace `out' = 3 if missing(`out') & `varlist'==0 & anti_hbs_pos==1 & anti_hbc_pos==1
    replace `out' = 4 if missing(`out') & `varlist'==0 & anti_hbs_pos==0 & anti_hbc_pos==1
    replace `out' = 5 if missing(`out') & `varlist'==0 & anti_hbs_pos==0 & anti_hbc_pos==0
    label values `out' piimm
end

* >1 categories (keeps your original but with explicit name)
__make_pi_cat hbsag_gt1,  out(pi_imm_cat_gt1)
label var pi_imm_cat_gt1 "PI schematic category (HBsAg S/CO>1)"

* >10 categories (new)
__make_pi_cat hbsag_gt10, out(pi_imm_cat_gt10)
label var pi_imm_cat_gt10 "PI schematic category (HBsAg S/CO>10)"

* (Optional sanity checks)
tab pi_imm_cat_gt1, m
tab pi_imm_cat_gt10, m

/******************************************************************
  E) Export tidy prevalence tables (overall) for BOTH thresholds
******************************************************************/
preserve
tempfile pi_prev_gt1 pi_prev_gt10
capture postclose P1
postfile P1 str50 category double prev se l u long N using `pi_prev_gt1', replace

forvalues k = 1/5 {
    quietly count if pi_imm_cat_gt1==`k'
    if (r(N)==0) continue
    tempvar z
    gen byte `z' = (pi_imm_cat_gt1==`k') if !missing(pi_imm_cat_gt1)
    quietly svy: mean `z'
    scalar b  = _b[`z']
    scalar se = _se[`z']
    scalar df = e(df_r)
    scalar t  = invttail(df, 0.025)
    scalar l  = b - abs(t)*se
    scalar u  = b + abs(t)*se
    quietly count if pi_imm_cat_gt1==`k'
    scalar N = r(N)
    local lab : label piimm `k'
    if "`lab'"=="" local lab "Category `k'"
    post P1 ("`lab'") (b) (se) (l) (u) (N)
}
postclose P1

use `pi_prev_gt1', clear
gen prev_pct = 100*prev
gen l_pct = 100*l
gen u_pct = 100*u
format prev_pct l_pct u_pct %6.1f
order category prev_pct l_pct u_pct prev se l u N
sort category
export excel using "$TABLES/pi_immunity_prevalence_gt1.xlsx", firstrow(variables) replace
display as text "Exported: $TABLES/pi_immunity_prevalence_gt1.xlsx"
restore
* ---- Now for >10
capture postclose P2
postfile P2 str50 category double prev se l u long N using `pi_prev_gt10', replace

forvalues k = 1/5 {
    quietly count if pi_imm_cat_gt10==`k'
    if (r(N)==0) continue
    tempvar z
    gen byte `z' = (pi_imm_cat_gt10==`k') if !missing(pi_imm_cat_gt10)
    quietly svy: mean `z'
    scalar b  = _b[`z']
    scalar se = _se[`z']
    scalar df = e(df_r)
    scalar t  = invttail(df, 0.025)
    scalar l  = b - abs(t)*se
    scalar u  = b + abs(t)*se
    quietly count if pi_imm_cat_gt10==`k'
    scalar N = r(N)
    local lab : label piimm `k'
    if "`lab'"=="" local lab "Category `k'"
    post P2 ("`lab'") (b) (se) (l) (u) (N)
}
postclose P2

use `pi_prev_gt10', clear
gen prev_pct = 100*prev
gen l_pct = 100*l
gen u_pct = 100*u
format prev_pct l_pct u_pct %6.1f
order category prev_pct l_pct u_pct prev se l u N
sort category
export excel using "$TABLES/pi_immunity_prevalence_gt10.xlsx", firstrow(variables) replace
display as text "Exported: $TABLES/pi_immunity_prevalence_gt10.xlsx"



*======================
* for the plot
*======================

* From ev_wd2.dta -> export the fields we need for plotting in R
use "$DERIVED/vk_evolve.dta", clear
svyset [pweight=totwt], singleunit(centered)

* Parse numeric anti-HBs
capture drop hbsab_num
gen double hbsab_num = .
replace hbsab_num = real(regexs(1)) if regexm(lower(hbsab_result),"([0-9.]+)")

* Seromarkers
gen byte hbsag_neg = lower(hbsag_interpretation)=="nonreactive"
gen byte anti_hbc_pos = lower(anti_hbcii_interpretation)=="reactive"

* Keep HBsAg-negative only (immunity context)
keep if hbsag_neg==1

* HIV labels (0=Negative,1=Positive assumed)
label define HIV 0 "HIV-negative" 1 "PLWH", replace
label values hivstat HIV

* Tidy export (drop missing or zero titres; R will handle log safely if ≥1)
keep iintid hbsag_neg anti_hbc_pos hivstat hbsab_num
drop if missing(hbsab_num)
replace hbsab_num = 1 if hbsab_num<1 & !missing(hbsab_num)   // floor for log-scale visibility
export delimited using "$OUT/immunity_titres_for_plot_hiv.csv", replace
di as res "Wrote: $OUT/immunity_titres_for_plot_hiv.csv"

use "$DERIVED/vk_evolve.dta", clear
svyset [pweight=totwt], singleunit(centered)

capture drop hbsab_num
gen double hbsab_num = .
replace hbsab_num = real(regexs(1)) if regexm(lower(hbsab_result),"([0-9.]+)")

gen byte hbsag_neg  = lower(hbsag_interpretation)=="nonreactive"
gen byte anti_hbc_pos = lower(anti_hbcii_interpretation)=="reactive"

keep if hbsag_neg==1
keep iintid anti_hbc_pos hbsab_num
drop if missing(hbsab_num)
replace hbsab_num = 1 if hbsab_num<1 & !missing(hbsab_num)

export delimited using "$OUT/immunity_titres_for_plot.csv", replace
