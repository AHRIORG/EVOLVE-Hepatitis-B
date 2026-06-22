/**********************************************************************
 EVOLVE HBV STUDY – Main analysis
 Author			: L. Mazibuko
 Date started 	: 26 February 2024
 Last updated	: 8 Sept 2025
Purpose:
   1) Link lab results to Vukuzazi core data and weights
   2) Construct outcomes:
        - HBsAg positive (manufacturer cut-off: S/CO > 1)
        - HBsAg strict sensitivity (S/CO >2, >5,  > 10, >100)
        - Vaccine-mediated immunity, Exposure & clearance, Susceptible
   3) Estimate weighted prevalence overall & stratified
   4) Fit prespecified multivariable models (no univariate fishing)
   5) Treatment eligibility: **APPLY ONLY to S/CO > 10** (confirmed)
***********************************************************************/

version 18
clear all
set more off
set linesize 255

/*==============================*
 | 0. Paths & logging
 *==============================*/
cd "/Users/lusandamazibuko/Library/CloudStorage/OneDrive-AHRI/Documents/EVOLVE Vukuzazi/evolve-hepatitis-b"
do "src/stata/00_preflight.do"
do "src/stata/00_run_all.do"
do "src/stata/_config.do"

capture log close _all
log using "$LOGS/evolve_hbv_main.smcl", replace

/*==============================*
 | 1. Input data
 *==============================*/

/* 1A. Core Vukuzazi dataset (from Glory)
   Keep in Stata format and standardise names
*/
use "$DATA/Vukuzazi_mortality_analysis.dta", clear
rename _all, lower
tempfile vukuzazi lab vk_main hbv_cases hbv_ctrl hbvflag
save `vukuzazi', replace

/* 1B. Lab results (latest processed Excel from March 2024)
   Import and standardise names
*/
import excel using "$DATA/VUKUZAZI EVOLVE HBV STUDY RESULTS_13MARCH2024_MOD.xlsx", ///
    sheet("Sheet1") firstrow clear
rename _all, lower
save `lab', replace

/*==============================*
 | 2. Merge lab + Vukuzazi
 *==============================*/

/* Ensure keys are unique before 1:1 merge */
isid individualid iintid using `lab', sort
isid individualid iintid using `vukuzazi', sort

use `lab', clear
merge 1:1 individualid iintid using `vukuzazi', keepusing(*) keep(match) nogen

/*==============================*
 | 3. Merge survey weights
 *==============================*/
merge 1:1 individualid using "$DERIVED/EVOLVE_weights.dta", keep(match) nogen

/*==============================*
 | 4. Variable hygiene & labels
 *==============================*/

/* Harmonise text for interpretations (guard against case/typos) */
foreach v in hbsag_interpretation hbsab_interpretation anti_hbcii_interpretation {
    replace `v' = strproper(strtrim(`v'))
}

/* Basic value labels */
capture label drop yesno
label define yesno 0 "No" 1 "Yes"

/* Sociodemographic labels */
label var sex           "Sex"
label var sescat        "Socioeconomic status"
label var drink3cat     "Alcohol intake"
label var smoke3cat     "Smoking status"
label var bmicat        "BMI category"

capture label drop smoke_ever_lbl
label define smoke_ever_lbl 0 "Never" 1 "Ever smoked (current/former)"

capture drop smoke_ever
gen byte smoke_ever = .
capture confirm variable smoke3cat
if !_rc {
    capture confirm string variable smoke3cat
    if !_rc {
        replace smoke_ever = 0 if regexm(lower(strtrim(smoke3cat)), "never")
        replace smoke_ever = 1 if regexm(lower(strtrim(smoke3cat)), "current|ex-smoker|ex smoker|former|ever")
    }
    else {
        replace smoke_ever = 1 if inlist(smoke3cat, 1, 2)
        replace smoke_ever = 0 if smoke3cat == 3
        count if missing(smoke_ever) & !missing(smoke3cat)
        if r(N) > 0 {
            tempvar smoke3cat_text
            capture decode smoke3cat, gen(`smoke3cat_text')
            if !_rc {
                replace smoke_ever = 0 if missing(smoke_ever) & regexm(lower(strtrim(`smoke3cat_text')), "never")
                replace smoke_ever = 1 if missing(smoke_ever) & regexm(lower(strtrim(`smoke3cat_text')), "current|ex-smoker|ex smoker|former|ever")
            }
        }
        count if missing(smoke_ever) & !missing(smoke3cat)
        if r(N) > 0 {
            di as error "WARNING: smoke3cat has unexpected values; some smoke_ever values remain missing."
        }
    }
}
else {
    di as error "WARNING: smoke3cat not found; smoke_ever remains missing."
}
label values smoke_ever smoke_ever_lbl
label var smoke_ever "Smoking: ever smoked (current/ex) vs never"

/* Education (ensure consistent numeric coding) */
capture confirm variable educ_socdem
if _rc==0 {
    encode educ_socdem, gen(educ_num)
    recode educ_num (1=0 "None") (2=1 "Primary") (3=2 "Secondary") (4=3 "Tertiary"), gen(educ_recoded)
    label define educ_lbl 0 "None" 1 "Primary" 2 "Secondary" 3 "Tertiary"
    label values educ_recoded educ_lbl
    label var educ_recoded "Education level"
}
else {
    di as err "NOTE: educ_socdem not found; using existing education variable(s) as-is."
}

/*==============================*
 | 5. Outcomes & key covariates
 *==============================*/

/* A) Primary outcomes based on qualitative interpretations
   Active infection: HBsAg reactive
   Vaccine-mediated immunity: HBsAb reactive, HBsAg nonreactive, anti-HBc (total) nonreactive
   Exposure & clearance: anti-HBc reactive, HBsAg nonreactive (HBsAb may be +/-)
   Susceptible: HBsAg-, anti-HBc-, HBsAb-
   
   
   anti-HBc / anti_hbcii_interpretation - a positive result for total anti-HBc, including from an Anti-HBc II test, indicates either a past or current Hepatitis B infection 
   and serves as a key marker to document exposure to the Hepatitis B virus (HBV). 
   
  anti-HBs OR HBsAb, or Hepatitis B surface antibody, is an antibody produced by the immune system that indicates protection against HBV. 
  A positive HBsAb result shows you are immune, either from recovering from a natural infection or from receiving a hepatitis B vaccine. 
  
  HBsAg is a protein found on the outer surface of the hepatitis B virus 
  
*/
/* Make named cleaned copies (kept with readable names in saved datasets) */
capture drop hbsag_clean hbsab_clean antihbc_clean
capture drop miss_hbsag miss_hbsab miss_ahbc

gen str30 hbsag_clean   = lower(trim(hbsag_interpretation))
gen str30 hbsab_clean   = lower(trim(hbsab_interpretation))
gen str30 antihbc_clean = lower(trim(anti_hbcii_interpretation))

/* A small helper to flag missing/insufficient/empty */
gen byte miss_hbsag = inlist(hbsag_clean,"","insufficient","equivocal")
gen byte miss_hbsab = inlist(hbsab_clean,"","insufficient","equivocal")
gen byte miss_ahbc  = inlist(antihbc_clean,"","insufficient","equivocal")

/* ====== Outcomes per protocol ====== */

/* (i) HBV infection (HBsAg positive, irrespective of other markers) */
capture drop active_disease
gen byte active_disease = .
replace active_disease = 1 if hbsag_clean=="reactive"
replace active_disease = 0 if hbsag_clean=="nonreactive"
label define yesno 0 "No" 1 "Yes", replace
label values active_disease yesno
label var     active_disease "HBsAg positive (qualitative)"

/* (ii) Exposure with clearance: HBsAg- & anti-HBc+ (anti-HBs can be +/-) */ ** before this was hbsab_interpretation negative instead of hbsag_interpretation
capture drop hbvexpo_clear
gen byte hbvexpo_clear = .
replace hbvexpo_clear = 1 if hbsag_clean=="nonreactive" & antihbc_clean=="reactive"
replace hbvexpo_clear = 0 if hbsag_clean=="nonreactive" & antihbc_clean=="nonreactive" ///
                            & !miss_ahbc
/* leave as . if any required marker missing/insufficient */
label values hbvexpo_clear yesno
label var     hbvexpo_clear "Exposure & clearance (anti-HBc+ / HBsAg-)"

/* (iii) Vaccine-mediated immunity: anti-HBs+ AND (HBsAg- & anti-HBc-) */
capture drop vax_med_immune
gen byte vax_med_immune = .
replace vax_med_immune = 1 if hbsab_clean=="reactive"  ///
                            & hbsag_clean=="nonreactive" ///
                            & antihbc_clean =="nonreactive"
replace vax_med_immune = 0 if hbsag_clean=="nonreactive" & antihbc_clean=="nonreactive" ///
                            & hbsab_clean=="nonreactive"
/* leave as . if any required marker missing/insufficient */
label values vax_med_immune yesno
label var     vax_med_immune "Vaccine-mediated immunity (HBsAb+ / others -)"

/* Susceptible: all three negative */
capture drop susceptible
gen byte susceptible = .
replace susceptible = 1 if hbsag_clean=="nonreactive" & antihbc_clean=="nonreactive" & hbsab_clean=="nonreactive"
replace susceptible = 0 if susceptible==. & (hbsag_clean!="" | antihbc_clean!="" | hbsab_clean!="")
label values susceptible yesno
label var     susceptible "Susceptible (HBsAg-/anti-HBc-/anti-HBs-)"


/* Mutually exclusive status among those with all three markers observed */
count if inlist(active_disease,. ,0,1) & inlist(hbvexpo_clear,. ,0,1) & inlist(vax_med_immune,. ,0,1)
tab active_disease hbvexpo_clear if vax_med_immune==1, m
tab vax_med_immune hbvexpo_clear if active_disease==1, m

/* Everyone classifiable when all three markers present? */
count if miss_hbsag==0 & miss_hbsab==0 & miss_ahbc==0
tab1 active_disease hbvexpo_clear vax_med_immune susceptible if miss_hbsag==0 & miss_hbsab==0 & miss_ahbc==0, m

/* B) Threshold-based sensitivity outcome (S/CO > 10)
   Extract numeric S/CO from hbsag_result string if needed
*/
capture confirm numeric variable hbsag_result
if _rc {
    // Extract first numeric token (allows 1.07, 10, 5101.18)
    gen double hbv_num = real(regexs(1)) if regexm(hbsag_result, "([0-9.]+)")
}
else {
    gen double hbv_num = hbsag_result
}

gen byte hbv_active = .
replace hbv_active = 0 if hbv_num <= 10 & !missing(hbv_num)
replace hbv_active = 1 if hbv_num >  10
label values hbv_active yesno
label var     hbv_active "HBsAg S/CO > 10 (sensitivity definition)"

/* C) Birth cohort (HBV vaccine epoch) */
capture confirm numeric variable dateofbirth
if _rc { // assume daily Stata date in dateofbirth_s
    // no-op; user can adapt if date variable differs
}
gen int birthyear = year(dateofbirth)
gen byte hbvdob = .
replace hbvdob = 0 if birthyear < 1995
replace hbvdob = 1 if inrange(birthyear,1995,1999)
replace hbvdob = 2 if inrange(birthyear,2000,2005)
label define hbvdob 0 "Born <1995" 1 "1995–1999" 2 "2000–2005", replace
label values hbvdob hbvdob
label var     hbvdob "HBV vaccine age epoch"

/* D) HIV status + care cascade */
gen byte hivstat = .
replace hivstat = 1 if hivelisa==1   // positive
replace hivstat = 0 if hivelisa==2   // negative
label define hivstat 0 "Negative" 1 "Positive", replace
label values hivstat hivstat
label var     hivstat "HIV status (ELISA)"

gen byte vl_suppressed = .
replace vl_suppressed = 1 if hivstat==1 & vl<50
replace vl_suppressed = 0 if hivstat==1 & vl>=50 & vl<.
label define supp 0 "Not suppressed" 1 "Suppressed", replace
label values vl_suppressed supp
label var     vl_suppressed "Viral load <50 copies/mL (PLWH only)"

/* E) Hypertension & Diabetes */
gen byte hypertension = .
replace hypertension = 1 if htncascade>0 & htncascade<.
replace hypertension = 0 if htncascade==0
label values hypertension yesno
label var     hypertension "Hypertension"

gen byte diabetic = .
replace diabetic = 1 if dmdiag==1
replace diabetic = 0 if dmdiag==2
label values diabetic yesno
label var     diabetic "Diabetes"

/* F) On ART (prefer cascade over self-report), among PLWH */
gen byte onart = .
replace onart = 0 if hivstat==1 & inrange(hivcascade,0,2)
replace onart = 1 if hivstat==1 & hivcascade>2 & hivcascade<.
label values onart yesno
label var     onart "On ART (cascade)"


/*============================================
Numeric value of VACCINE-MEDIATED IMMUNITY
==============================================
*/
*--- Clean the raw string and normalize units
capture drop hbsab_result_clean
gen strL hbsab_result_clean = lower(trim(hbsab_result))
replace hbsab_result_clean = subinstr(hbsab_result_clean,"miu/ml","",.)
replace hbsab_result_clean = subinstr(hbsab_result_clean,"mIU/mL","",.)   // just in case mixed case
replace hbsab_result_clean = strtrim(hbsab_result_clean)

*--- Capture censoring operators (<, >) and extract the numeric part
gen byte hbsab_op = .              // -1="<", 0="exact", +1=">"
replace hbsab_op = -1 if regexm(hbsab_result_clean,"^\s*<")
replace hbsab_op =  1 if regexm(hbsab_result_clean,"^\s*>")
replace hbsab_op =  0 if regexm(hbsab_result_clean,"^[0-9]")

gen double hbsab_num = .
replace hbsab_num = real(regexs(1)) if regexm(hbsab_result_clean,"([0-9]+(\.[0-9]+)?)")

label define HBSAB_OP -1 "< (left-censored)" 0 "exact" 1 "> (right-censored)"
label values hbsab_op HBSAB_OP
label var hbsab_num "anti-HBs (mIU/mL), numeric (raw parsed)"
label var hbsab_op  "anti-HBs censoring flag"

*--- Protection indicator at 10 mIU/mL
gen byte hbsab_ge10 = (hbsab_num >= 10) if hbsab_num < .
label define YESNO 0 "No" 1 "Yes", replace
label values hbsab_ge10 YESNO
label var hbsab_ge10 "anti-HBs ≥10 mIU/mL"

*======
* PLOT 
*======
* keep only nonmissing numeric values
keep if !missing(hbsab_num)

* x = observation index
gen long obs = _n

twoway scatter hbsab_num obs, ///
    msize(vsmall) msymbol(o) ///
    ytitle("anti-HBs (mIU/mL)") xtitle("Observation") ///
    title("anti-HBs values by observation") legend(off)


twoway scatter hbsab_num obs, ///
    msize(vsmall) msymbol(o) yscale(log) ///
    ytitle("anti-HBs (mIU/mL, log scale)") xtitle("Observation") ///
    title("anti-HBs (log scale)") legend(off)

	
	
	
	
/*==============================*
 | 6. Save analysis dataset
 *==============================*/
order individualid iintid sex birthyear hbvdob ///
      educ_recoded sescat drink3cat smoke_ever smoke3cat bmicat ///
      hivstat onart vl_suppressed ///
      hypertension diabetic ///
      hbsag_interpretation hbsab_interpretation anti_hbcii_interpretation ///
      hbsag_result hbv_num active_disease hbv_active vax_med_immune hbvexpo_clear susceptible ///
      ipw_weights totwt

drop if missing(totwt)
drop_stata_tempvars
save "$DERIVED/vk_evolve.dta", replace

/*==============================*
 | 7. Survey design
 *==============================*/
use "$DERIVED/vk_evolve.dta", clear
svyset [pweight=totwt]

/*==============================*
 | 8. Table 1 – weighted descriptives
 *==============================*

cap which summtab
if _rc==0 {
    summtab, contvars(ageatenrolment) ///
             catvars(hbvdob sex educ_recoded sescat drink3cat smoke_ever ///
                     hivstat onart vl_suppressed bmicat hypertension diabetic) ///
             mean median total ///
             title("Table 1. Sociodemographic & clinical characteristics (weighted)") ///
             excel excelname("$OUT/summary_table.xlsx") replace
}
else {
    // Fallback: quick weighted one-way proportions to Excel
    postutil clear
    tempfile t1
    postfile H str32 var level double p se using `t1', replace
    foreach v of varlist hbvdob sex educ_recoded sescat drink3cat smoke_ever ///
                         hivstat onart vl_suppressed bmicat hypertension diabetic {
        levelsof `v', local(L)
        foreach L of local L {
            quietly svy: proportion i.`v' if `v'==`L'
            matrix M = r(table)
            local p  = M[1,1]
            local se = M[2,1]
            post H ("`v'") ("`L'") (`p') (`se')
        }
    }
    postclose H
    use `t1', clear
    export excel using "$OUT/summary_table_fallback.xlsx", firstrow(variables) replace
}
*/

/*==============================*
 | 9. Table 2 – weighted prevalence (overall & stratified)
 *==============================*/
svyset [pweight=totwt]

/* Outcomes: active_disease, hbv_active (S/CO≥10), vax_med_immune, hbvexpo_clear, susceptible */
svy: proportion active_disease vax_med_immune hbvexpo_clear susceptible hbv_active

/* Stratified by sex, birth cohort, HIV status */
foreach s in sex hbvdob hivstat {
    di as txt "---- Weighted prevalence by `s' ----"
    svy, over(`s'): proportion active_disease vax_med_immune hbvexpo_clear susceptible hbv_active
}

/* Note: `svy: tabulate` supports row/column % and SEs; Pearson chi2 is not appropriate under svy.
         If you want tests, use `svy: tabulate var outcome, row se` and interpret design-based F. */
svyset [pweight=totwt]

* --- Overall (no p-values) ---
svy: proportion active_disease vax_med_immune hbvexpo_clear

* --- Design-based p-values by SEX ---
svy: tabulate sex active_disease, row pearson     // reports design-based F and p
svy: tabulate sex vax_med_immune, row pearson
svy: tabulate sex hbvexpo_clear, row pearson

* --- Design-based p-values by HBV vaccine age epoch (hbvdob) ---
svy: tabulate hbvdob active_disease, row pearson
svy: tabulate hbvdob vax_med_immune, row pearson
svy: tabulate hbvdob hbvexpo_clear, row pearson

* --- Design-based p-values by HIV status ---
svy: tabulate hivstat active_disease, row pearson
svy: tabulate hivstat vax_med_immune, row pearson
svy: tabulate hivstat hbvexpo_clear, row pearson

		 
		 
/*===============================================================*
 | Threshold grid (strict ">"): >1, >2, >5, >10, >100
 | Outputs: weighted prevalence overall and by hbvdob, sex, hivstat
 | Exports: Excel at $OUT/hbsag_threshold_prevalence.xlsx
 *===============================================================*/

use "$DERIVED/vk_evolve.dta", clear

/* Survey design */
svyset [pweight=totwt]
svyset, singleunit(centered)   // keeps CIs defined in tiny subpops

/* Ensure hbv_num exists (numeric S/CO) */
capture confirm numeric variable hbv_num
if _rc {
    gen double hbv_num = real(regexs(1)) if regexm(hbsag_result, "([0-9.]+)")
}

/* 1) STRICT indicators: POS = (S/CO > k), NEG = (S/CO ≤ k) */
local cuts 1 2 5 10 100
label define YESNO 0 "No" 1 "Yes", replace

foreach k of local cuts {
    /* generate (or replace safely) */
    capture drop hbsag_gt`k' hbsag_le`k'
    gen byte hbsag_gt`k' = (hbv_num >  `k') if !missing(hbv_num)   // POSITIVE
    gen byte hbsag_le`k' = (hbv_num <= `k') if !missing(hbv_num)   // NEGATIVE
    label values hbsag_gt`k' YESNO
    label values hbsag_le`k' YESNO
    label var hbsag_gt`k' "HBsAg S/CO > `k' (POS)"
    label var hbsag_le`k' "HBsAg S/CO <= `k' (NEG)"
}

/* --- DEFINE helper with a proper subpop flag --- */
cap program drop __post_prev
program define __post_prev, rclass
    version 18
    syntax varname [if], THRESH(string) SIDE(string) GROUPVAR(string) LEVEL(string)

    tempvar sp
    /* If no IF was passed (Overall), set subpop=1 for all; otherwise flag the subset */
    if "`if'"=="" {
        gen byte `sp' = 1
    }
    else {
        gen byte `sp' = 0
        quietly replace `sp' = 1 `if'
    }

    quietly svy, subpop(`sp'): mean `varlist'
    local b  = _b[`varlist']
    local se = _se[`varlist']
    local df = e(df_r)
    local t  = invttail(`df', 0.025)
    local l  = `b' - abs(`t')*`se'
    local u  = `b' + abs(`t')*`se'

    /* Unweighted N in the subpopulation (for context) */
    quietly count if `sp'==1
    local n = r(N)

    /* Store for potential p-value calculation */
    return local b = `b'
    return local se = `se'
    
    post P ("`thresh'") ("`side'") ("`groupvar'") ("`level'") (`b') (`se') (`l') (`u') (`n') (.)
end

/* --- Helper for group comparisons with p-values --- */
cap program drop __post_prev_comp
program define __post_prev_comp
    version 18
    syntax varname [if], THRESH(string) SIDE(string) GROUPVAR(string) LEVEL(string) REFCOND(string)

    tempvar sp sp_ref
    
    /* Create subpop indicators */
    if "`if'"=="" {
        gen byte `sp' = 1
    }
    else {
        gen byte `sp' = 0
        quietly replace `sp' = 1 `if'
    }
    
    gen byte `sp_ref' = 0
    quietly replace `sp_ref' = 1 `refcond'
    
    /* Get estimates for current group */
    quietly svy, subpop(`sp'): mean `varlist'
    local b  = _b[`varlist']
    local se = _se[`varlist']
    local df = e(df_r)
    local t  = invttail(`df', 0.025)
    local l  = `b' - abs(`t')*`se'
    local u  = `b' + abs(`t')*`se'
    
    quietly count if `sp'==1
    local n = r(N)
    
    /* Get reference group estimate */
    quietly svy, subpop(`sp_ref'): mean `varlist'
    local b_ref = _b[`varlist']
    local se_ref = _se[`varlist']
    
    /* Calculate p-value using z-test for difference in proportions */
    local diff = `b' - `b_ref'
    local se_diff = sqrt(`se'^2 + `se_ref'^2)
    local z = `diff' / `se_diff'
    local pval = 2 * (1 - normal(abs(`z')))
    
    post P ("`thresh'") ("`side'") ("`groupvar'") ("`level'") (`b') (`se') (`l') (`u') (`n') (`pval')
end

/* === BUILD THE RESULTS TABLE === */
tempfile T
postutil clear
postfile P str8 threshold str5 side str12 groupvar str20 level ///
        double(prev se ci_l ci_u n_eff pval) using `T', replace

/* 3) Overall + stratified, BOTH sides (>k POS, ≤k NEG) */
local cuts 1 2 5 10 100
foreach k of local cuts {
    foreach side in gt le {
        local var hbsag_`side'`k'
        
        /* Overall (no p-value) */
        __post_prev `var', thresh(">`k'") side("`=cond("`side'"=="gt","POS","NEG")'") ///
            groupvar("Overall") level("All")
        
        /* hbvdob: compare to reference (Born <1995 = 0) */
        levelsof hbvdob, local(Lhbv)
        foreach L of local Lhbv {
            if `L' == 0 {
                __post_prev `var' if hbvdob==`L', thresh(">`k'") side("`=cond("`side'"=="gt","POS","NEG")'") ///
                    groupvar("hbvdob") level("`L'")
            }
            else {
                __post_prev_comp `var' if hbvdob==`L', thresh(">`k'") side("`=cond("`side'"=="gt","POS","NEG")'") ///
                    groupvar("hbvdob") level("`L'") refcond("if hbvdob==0")
            }
        }
        
        /* sex: compare to reference (Female = 2) */
        levelsof sex, local(Lsex)
        foreach L of local Lsex {
            if `L' == 2 {
                __post_prev `var' if sex==`L', thresh(">`k'") side("`=cond("`side'"=="gt","POS","NEG")'") ///
                    groupvar("sex") level("`L'")
            }
            else {
                __post_prev_comp `var' if sex==`L', thresh(">`k'") side("`=cond("`side'"=="gt","POS","NEG")'") ///
                    groupvar("sex") level("`L'") refcond("if sex==2")
            }
        }
        
        /* hivstat: compare to reference (HIV-negative = 0) */
        levelsof hivstat, local(Lhiv)
        foreach L of local Lhiv {
            if `L' == 0 {
                __post_prev `var' if hivstat==`L', thresh(">`k'") side("`=cond("`side'"=="gt","POS","NEG")'") ///
                    groupvar("hivstat") level("`L'")
            }
            else {
                __post_prev_comp `var' if hivstat==`L', thresh(">`k'") side("`=cond("`side'"=="gt","POS","NEG")'") ///
                    groupvar("hivstat") level("`L'") refcond("if hivstat==0")
            }
        }
    }
}

/* 4) Finalise tidy table */
postclose P
use `T', clear
format prev se ci_l ci_u %9.4f
format pval %6.4f

gen level_label = level
replace level_label = "Born <1995"   if groupvar=="hbvdob" & level=="0"
replace level_label = "1995–1999"    if groupvar=="hbvdob" & level=="1"
replace level_label = "2000–2005"    if groupvar=="hbvdob" & level=="2"
replace level_label = "Female"       if groupvar=="sex"     & level=="2"
replace level_label = "Male"         if groupvar=="sex"     & level=="1"
replace level_label = "HIV-negative" if groupvar=="hivstat" & level=="0"
replace level_label = "HIV-positive" if groupvar=="hivstat" & level=="1"
replace level_label = "All"          if groupvar=="Overall"

gen prev_pct = 100*prev
gen l_pct    = 100*ci_l
gen u_pct    = 100*ci_u
format prev_pct l_pct u_pct %6.2f

/* Add significance stars */
gen sig_star = ""
replace sig_star = "***" if pval < 0.001 & !missing(pval)
replace sig_star = "**"  if pval < 0.01  & pval >= 0.001 & !missing(pval)
replace sig_star = "*"   if pval < 0.05  & pval >= 0.01  & !missing(pval)

order threshold side groupvar level_label prev_pct l_pct u_pct pval sig_star n_eff
sort threshold side groupvar level_label

drop_stata_tempvars
save "$DERIVED/hbsag_threshold_prevalence_long.dta", replace

/* 5) Export tidy long sheet */
drop_stata_tempvars
export excel using "$OUT/hbsag_threshold_prevalence.xlsx", firstrow(variables) replace

/* 6) Multi-sheet workbook */
local book "$OUT/hbsag_threshold_prevalence_bygroup.xlsx"
capture erase "`book'"

use "$DERIVED/hbsag_threshold_prevalence_long.dta", clear
foreach g in Overall hbvdob sex hivstat {
    preserve
    keep if groupvar=="`g'"
    if "`g'"=="Overall" {
        export excel using "`book'", sheet("`g'", replace) firstrow(variables)
    }
    else {
        export excel using "`book'", sheet("`g'", modify) firstrow(variables)
    }
    restore
}

display as result _n "Analysis complete!"
display as result "Results saved to:"
display as text "  - $OUT/hbsag_threshold_prevalence_long.dta"
display as text "  - $OUT/hbsag_threshold_prevalence.xlsx"
display as text "  - $OUT/hbsag_threshold_prevalence_bygroup.xlsx"

/*==============================*
 | 9b. Table 3 – prespecified multivariable models (weighted)
 *==============================*/
use "$DERIVED/vk_evolve.dta", clear

/* One prespecified model per outcome (no separate "univariate" regressions) */
global ADJ i.hbvdob i.sex i.educ_recoded i.sescat i.drink3cat i.smoke_ever ///
          i.hivstat ib2.bmicat i.hypertension i.diabetic

eststo clear
svyset [pweight=totwt]


svy: logistic active_disease $ADJ
eststo m_active
svy: logistic hbv_active    $ADJ
eststo m_active10
svy: logistic vax_med_immune $ADJ
eststo m_vax
svy: logistic hbvexpo_clear  $ADJ
eststo m_clear

/* Export ORs with 95% CI 
esttab m_active m_active10 m_vax m_clear using "$OUT/table3_models.rtf", ///
    eform b(3) ci(3) p(3) nogap onecell star(* 0.05 ** 0.01 *** 0.001) ///
    title("Table 3. Adjusted odds ratios (survey-weighted)") replace
*/ 
/* Optional --- marginal prevalence by key factors for figures */
margins, at(hbvdob=(0 1 2)) post
marginsplot, name(mp_hbvdob, replace)

/*==============================*
 | 10. Case subset data for treatement eligibility assessment 
 *==============================*/

/* Cases (lab-confirmed set file) */
* NOTE: We will still restrict treatment eligibility to S/CO > 10 later.
use "$DERIVED/vk_evolve.dta", clear
drop_stata_tempvars
save `vk_main', replace 

import excel using "$DATA/EvoLVE AHRI_PLA_for HBV_VL-HBeAg-Alt_RESULTS 18 Oct 2024 1 Mod.xlsx", ///
    sheet("Sheet1") firstrow clear
rename _all, lower
save `hbv_cases', replace 

use `vk_main', clear
merge 1:1 iintid using `hbv_cases', nogen keep(match)
gen byte case = 1
label values case YESNO
drop_stata_tempvars
save "$DERIVED/hbv_cases.dta", replace

/* Controls */
import excel using "$DATA/AHRI_SER_for Alt_RESULTS 2025-01-25 Mod.xlsx", ///
    sheet("Sheet1") firstrow clear
rename _all, lower
save `hbv_ctrl', replace 

use `hbv_ctrl', clear
merge 1:1 iintid using `vk_main', nogen keep(match)
gen byte case = 0
label values case YESNO
drop_stata_tempvars
save "$DERIVED/hbv_controls.dta", replace

/* Append cases + controls */
use "$DERIVED/hbv_controls.dta", clear
append using "$DERIVED/hbv_cases.dta"
label define CS 0 "Control" 1 "Case", replace
gen byte status = case
label values status CS
label var status "Case or control status"

svyset [pweight=totwt]

/*==============================*
 | 11. Virology & ALT handling (treatment eligibility)
 *==============================*/

/* HBV DNA VL harmonisation */
gen str20 dna_vl_str = hbv_viral_load_iu_ml
replace dna_vl_str = "9.9"         if hbv_viral_load_iu_ml=="<10"
replace dna_vl_str = "1000000001"  if hbv_viral_load_iu_ml==">1000000000"
replace dna_vl_str = "0"           if hbv_viral_load_iu_ml=="Target Not Detected"

destring dna_vl_str, gen(dna_vl) force
label var dna_vl "HBV DNA IU/mL (imputed bounds for <10, >1e9, TND)"

gen byte dna_thres = .
replace dna_thres = 1 if dna_vl < 10
replace dna_thres = 0 if dna_vl >=10 & dna_vl<.
label define DNACAT 0 "Not suppressed (>=10)" 1 "Suppressed (<10)", replace
label values dna_thres DNACAT

/* log-scale summaries (geometric means) – optional */
gen double logvl = .
replace logvl = ln(dna_vl) if dna_vl>0

/* ALT harmonisation */
replace alt_u_l = subinstr(alt_u_l,"<7,00","6.99",.)
replace alt_u_l = subinstr(alt_u_l,"<7","6.99",.)
destring alt_u_l, replace force
label var alt_u_l "ALT (U/L)"

/* WHO ULN categories (primary) */
gen byte ALT_category = .
replace ALT_category = 1 if (sex==1 & alt_u_l>30) | (sex==2 & alt_u_l>19) // Elevated
replace ALT_category = 0 if missing(ALT_category) & alt_u_l<.            // Normal
label define ALTCAT 0 "Normal (WHO)" 1 "Elevated (WHO)", replace
label values ALT_category ALTCAT
label var ALT_category "ALT category (WHO ULN)"

/* NHLS ULN categories (supplementary) */
gen byte nhls_altcat = .
replace nhls_altcat = 1 if (sex==1 & alt_u_l>40) | (sex==2 & alt_u_l>35) // Elevated
replace nhls_altcat = 0 if missing(nhls_altcat) & alt_u_l<.
label define ALTNHLS 0 "Normal (NHLS)" 1 "Elevated (NHLS)", replace
label values nhls_altcat ALTNHLS
label var nhls_altcat "ALT category (NHLS ULN)"

/* HBeAg status */
gen byte HBeAg_status = .
replace HBeAg_status = 1 if lower(hbeag_intepretation)=="reactive"
replace HBeAg_status = 0 if lower(hbeag_intepretation)=="nonreactive"
label define HBE 0 "HBeAg Negative" 1 "HBeAg Positive", replace
label values HBeAg_status HBE
label var HBeAg_status "HBeAg status"

svy: proportion ALT_category nhls_altcat HBeAg_status if status==1
svy: proportion ALT_category nhls_altcat HBeAg_status

/*==============================*
 | 12. **Treatment eligibility**
 | Apply ONLY to confirmed infection: S/CO > 10
 *==============================*/

/* Confirmed infection flag from main data (hbv_active) */
* merge hbv_active if not already present (ensures consistency)
preserve
use "$DERIVED/vk_evolve.dta", clear
keep iintid hbv_active
save `hbvflag', replace 
restore
merge 1:1 iintid using `hbvflag', nogen keep(master match)

/* Strictly limit eligibility evaluation to confirmed infection */
gen byte treateligible = .                                // missing by default
replace treateligible = 0 if status==1 & hbv_active==1    // initialise among confirmed cases
replace treateligible = 1 if status==1 & hbv_active==1 & ///
    ( dna_vl>2000 | ALT_category==1 | diabetic==1 | (hivstat==1 & onart==0) | hypertension==1 )

label values treateligible YESNO
label var     treateligible "Eligible for treatment (WHO-like; applied ONLY if S/CO > 10)"
svy: proportion treateligible if status==1 & hbv_active==1

/* DNA categories (descriptive, among cases) */
gen byte hbvlcat = .
replace hbvlcat = 1 if status==1 & dna_vl<10
replace hbvlcat = 2 if status==1 & inrange(dna_vl,10,2000)
replace hbvlcat = 3 if status==1 & inrange(dna_vl,2001,20000)
replace hbvlcat = 4 if status==1 & dna_vl>20000
label define VLBIN 1 "<10" 2 "10–2000" 3 "2001–20000" 4 ">20000", replace
label values hbvlcat VLBIN
svy: proportion hbvlcat if status==1

/* Descriptives by ART status among cases */
svy, over(onart): proportion hbvlcat HBeAg_status if status==1
svy, over(onart): mean alt_u_l if status==1

/* Descriptives by HIV status among cases */
svy: proportion hbvlcat if status==1 & hivstat==0
svy, over(onart): proportion hbvlcat if status==1 & hivstat==1
svy: mean alt_u_l if status==1 & hivstat==0

/

/*==============================*
 | 13. Close
 *==============================*/
log close
display as text "Analysis complete. Outputs in: $OUT"
