
 /**********************************************************************
 EVOLVE HBV STUDY – Paper model table plus S/CO>=10 sensitivity
 Author         : Lusanda Sanda Mazibuko
 Purpose:
   1) Construct S/CO-based outcomes using HBsAg numeric S/CO
      - HBV infection (S/CO>1)
      - HBV infection (S/CO>10)  
      - Paper table outcomes for vaccine-mediated immunity, exposure & clearance,
        and susceptible use the confirmed-negative framework: HBsAg S/CO <10.
   2) Fit prespecified multivariable models (svy: logistic)
   3) Export a publication-style table with Reference rows

 ASSUMES in memory:
   - totwt
   - hbv_num (HBsAg S/CO numeric) OR hbsag_result to parse
   - hbsag_interpretation hbsab_interpretation anti_hbcii_interpretation
   - covariates:
       ageatenrolment hbvdob sex educ_recoded sescat drink3cat smoke_ever
       hivstat bmicat hypertension diabetic
**********************************************************************/

version 18
clear all
set more off
set linesize 255

*-------------------------------*
* Paths
*-------------------------------*
do "src/stata/_config.do"
cap mkdir "$OUT"

*-------------------------------*
* Load analysis dataset
*-------------------------------*
use "$DERIVED/vk_evolve.dta", clear

/**********************************************************************
 EVOLVE – paper table outcomes
**********************************************************************/

svyset [pweight=totwt]

capture label drop yesno
label define yesno 0 "No" 1 "Yes", replace

/**********************************************************************
 1) Numeric HBsAg S/CO
**********************************************************************/
capture confirm var hbv_num
if _rc {
    capture confirm numeric variable hbsag_result
    if _rc {
        gen double hbv_num = real(regexs(1)) if regexm(hbsag_result, "([0-9.]+)")
    }
    else {
        gen double hbv_num = hbsag_result
    }
}
label var hbv_num "HBsAg S/CO (numeric)"

/**********************************************************************
 2) Clean interpretation strings
**********************************************************************/
capture drop _hbsag_i _hbsab_i _ahbc_i _miss_hbsag _miss_hbsab _miss_ahbc
gen str20 _hbsag_i = lower(trim(hbsag_interpretation))
gen str20 _hbsab_i = lower(trim(hbsab_interpretation))
gen str20 _ahbc_i  = lower(trim(anti_hbcii_interpretation))

gen byte _miss_hbsag = inlist(_hbsag_i,"","insufficient","equivocal")
gen byte _miss_hbsab = inlist(_hbsab_i,"","insufficient","equivocal")
gen byte _miss_ahbc  = inlist(_ahbc_i ,"","insufficient","equivocal")

/**********************************************************************
 3) Outcomes
**********************************************************************/
capture drop hbv_inf_gt1 hbv_inf_ge10 hbv_neg_lt10
gen byte hbv_inf_gt1  = (hbv_num >  1) if hbv_num < .
gen byte hbv_inf_ge10 = (hbv_num > 10) if hbv_num < .
gen byte hbv_neg_lt10 = (hbv_num <=  10) if hbv_num < .
label values hbv_inf_gt1 hbv_inf_ge10 hbv_neg_lt10 yesno

capture drop vax_med_immune_ge10
gen byte vax_med_immune_ge10 = .
replace vax_med_immune_ge10 = 1 if hbv_neg_lt10==1 & _hbsab_i=="reactive" & _ahbc_i=="nonreactive"
replace vax_med_immune_ge10 = 0 if hbv_neg_lt10<. & !_miss_hbsab & !_miss_ahbc & missing(vax_med_immune_ge10)
label values vax_med_immune_ge10 yesno

capture drop hbvexpo_clear_ge10
gen byte hbvexpo_clear_ge10 = .
replace hbvexpo_clear_ge10 = 1 if hbv_neg_lt10==1 & _ahbc_i=="reactive"
replace hbvexpo_clear_ge10 = 0 if hbv_neg_lt10<. & !_miss_ahbc & missing(hbvexpo_clear_ge10)
label values hbvexpo_clear_ge10 yesno

capture drop susceptible_ge10
gen byte susceptible_ge10 = .
replace susceptible_ge10 = 1 if hbv_neg_lt10==1 & _ahbc_i=="nonreactive" & _hbsab_i=="nonreactive"
replace susceptible_ge10 = 0 if hbv_neg_lt10<. & !_miss_ahbc & !_miss_hbsab & missing(susceptible_ge10)
label values susceptible_ge10 yesno

capture confirm variable vax_med_immune
if _rc {
    gen byte vax_med_immune = .
    replace vax_med_immune = 1 if _hbsag_i=="nonreactive" & _ahbc_i=="nonreactive" & _hbsab_i=="reactive"
    replace vax_med_immune = 0 if !_miss_hbsag & !_miss_ahbc & !_miss_hbsab & missing(vax_med_immune)
    label values vax_med_immune yesno
}

capture confirm variable hbvexpo_clear
if _rc {
    gen byte hbvexpo_clear = .
    replace hbvexpo_clear = 1 if _hbsag_i=="nonreactive" & _ahbc_i=="reactive"
    replace hbvexpo_clear = 0 if !_miss_hbsag & !_miss_ahbc & missing(hbvexpo_clear)
    label values hbvexpo_clear yesno
}

capture confirm variable susceptible
if _rc {
    gen byte susceptible = .
    replace susceptible = 1 if _hbsag_i=="nonreactive" & _ahbc_i=="nonreactive" & _hbsab_i=="nonreactive"
    replace susceptible = 0 if !_miss_hbsag & !_miss_ahbc & !_miss_hbsab & missing(susceptible)
    label values susceptible yesno
}

capture confirm variable smoke_ever
if _rc {
    capture label drop smoke_ever_lbl
    label define smoke_ever_lbl 0 "Never" 1 "Ever smoked (current/former)"

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
}

/**********************************************************************
 4) Run and display each model with base levels shown
**********************************************************************/
local covars ""

foreach spec in ///
    "hbvdob:i." ///
    "sex:i." ///
    "educ_recoded:i." ///
    "sescat:i." ///
    "drink3cat:i." ///
    "hivstat:i." ///
    "bmicat:ib2." ///
    "hypertension:i." ///
    "diabetic:i." {
    gettoken v prefix : spec, parse(":")
    local prefix = substr("`prefix'", 2, .)
    capture confirm variable `v'
    if !_rc {
        local covars "`covars' `prefix'`v'"
    }
    else {
        di as text "NOTE: optional model covariate `v' not found; skipping."
    }
}

capture confirm variable smoke_ever
if !_rc {
    local covars "`covars' i.smoke_ever"
}
else {
    capture confirm variable smoke3cat
    if !_rc {
        di as text "NOTE: smoke_ever not found; falling back to i.smoke3cat."
        local covars "`covars' i.smoke3cat"
    }
    else {
        di as text "NOTE: optional model covariate smoke_ever/smoke3cat not found; skipping smoking."
    }
}

di as text "Model covariates: `covars'"

di as text _newline(2) "============================================================"
di as text "MODEL 1: HBV infection (HBsAg S/CO > 1)"
di as text "============================================================"
svy: logistic hbv_inf_gt1 `covars', baselevels

di as text _newline(2) "============================================================"
di as text "MODEL 2: HBV infection (HBsAg S/CO > 10)"
di as text "============================================================"
svy: logistic hbv_inf_ge10 `covars', baselevels

di as text _newline(2) "============================================================"
di as text "MODEL 3: Vaccine-mediated immunity (paper table; HBsAg S/CO < 10 framework)"
di as text "============================================================"
svy: logistic vax_med_immune_ge10 `covars', baselevels

di as text _newline(2) "============================================================"
di as text "MODEL 4: HBV exposure and clearance (paper table; HBsAg S/CO < 10 framework)"
di as text "============================================================"
svy: logistic hbvexpo_clear_ge10 `covars', baselevels

di as text _newline(2) "============================================================"
di as text "MODEL 5: Susceptible (paper table; HBsAg S/CO < 10 framework)"
di as text "============================================================"
svy: logistic susceptible_ge10 `covars', baselevels
