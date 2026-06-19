/********************************************************************
 EVOLVE HBV – Canonical dataset builder
 Author: Lusanda Mazibuko
 Date:   2025-10-05
 Stata:  18

 Outputs:
   1) $DATA/vk_master.dta              -- Vukuzazi + serology + weights (clean)
   2) $DATA/hbv_cc_merged.dta          -- Cases + controls appended with VK covars
   3) $DATA/analysis_eligibility.dta   -- Canonical, analysis-ready (PRIMARY)

 Notes:
   - Primary "confirmed infection" definition = HBsAg S/CO > 10 (strict >).
   - Sensitivity flags for >2 and >5 included (per PI request).
   - Treatment eligibility only evaluated among confirmed infection groups.
 
 Panels: [A] S/CO>1 with DNA+ALT+comorbidity, [B] S/CO>10 with DNA+ALT+comorbidity
         [C] S/CO>1 ALT+comorbidity only,   [D] S/CO>10 ALT+comorbidity only
 Key PI rules implemented:
   • Evaluate HBsAg positives using strict ">k" (NEG is ≤k).
   • "Already on treatment" = PLWH who are on ART (onart==1 & hivstat==1).
   • For A/B: first restrict to DNA-detectable (dna_vl>0) to avoid false positives,
              then eligible if (DNA>2000 & ALT>ULN) OR (PLWH not on ART) OR diabetes.
   • For C/D: ignore DNA entirely; eligible if ALT>ULN OR (PLWH not on ART) OR diabetes.
   • ALT ULN (WHO): >30 U/L for males, >19 U/L for females.
********************************************************************/

version 18
quietly clear all
set more off

*-------------------------------*
* 0. Paths
*-------------------------------*
do "src/stata/_config.do"
cap mkdir "$OUT"

*-------------------------------*
* 0b. Label scaffolding (define once)
*-------------------------------*
label define YESNO 0 "No" 1 "Yes", replace
label define HIV   0 "Negative" 1 "Positive", replace
label define HBD   0 "Born <1995" 1 "1995–1999" 2 "2000–2005", replace
label define ALTCAT   0 "Normal (WHO)" 1 "Elevated (WHO)", replace
label define ALTNHLS  0 "Normal (NHLS)" 1 "Elevated (NHLS)", replace
label define CS    0 "Control" 1 "Case", replace
label define VLBIN 1 "<10" 2 "10–2000" 3 "2001–20000" 4 ">20000", replace
label define SUPP  0 "Not suppressed" 1 "Suppressed", replace

*-------------------------------*
* 1. Vukuzazi + Serology + Weights -> vk_master.dta
*-------------------------------*

* 1.1 Vukuzazi core
use "$DATA/Vukuzazi_mortality_analysis.dta", clear
rename _all, lower
compress
tempfile vk
save `vk'

* 1.2 Lab serology (HBsAg/HBsAb/anti-HBc)
import excel using "$DATA/VUKUZAZI EVOLVE HBV STUDY RESULTS_13MARCH2024_MOD.xlsx", ///
    sheet("Sheet1") firstrow clear
rename _all, lower
compress
tempfile lab
save `lab'

* 1.3 Merge VK + Lab on IDs
use `lab', clear
cap confirm var iintid
if _rc {
    di as err "ERROR: iintid not found in lab sheet."
    exit 198
}
merge 1:1 individualid iintid using `vk', keep(match) nogen

* 1.4 Merge survey weights
merge 1:1 individualid using "$DERIVED/EVOLVE_weights.dta", keep(match) nogen
label var totwt "Survey weight (selection+nonresponse)"

*-------------------------------*
* 1.5 Core outcomes & covariates (single definitions)
*-------------------------------*

* Parse numeric S/CO (hbv_num) from HBsAg result (string or numeric)
cap confirm numeric variable hbsag_result
if _rc {
    gen double hbv_num = real(regexs(1)) if regexm(hbsag_result, "([0-9.]+)")
}
else {
    gen double hbv_num = hbsag_result
}
label var hbv_num "HBsAg S/CO (numeric)"

* Manufacturer qualitative HBsAg
capture drop active_disease
gen byte active_disease = .
replace active_disease = 1 if lower(hbsag_interpretation)=="reactive"
replace active_disease = 0 if lower(hbsag_interpretation)=="nonreactive"
label values active_disease YESNO
label var     active_disease "HBsAg positive (qualitative)"

* Derived serostatus
capture drop vax_med_immune
gen byte vax_med_immune = .
replace vax_med_immune = 1 if lower(hbsab_interpretation)=="reactive"  ///
                        & lower(hbsag_interpretation)=="nonreactive"   ///
                        & lower(anti_hbcii_interpretation)=="nonreactive"
replace vax_med_immune = 0 if vax_med_immune==.
label values vax_med_immune YESNO
label var     vax_med_immune "Vaccine-mediated immunity (HBsAb+ / HBsAg- / anti-HBc-)"

capture drop hbvexpo_clear
gen byte hbvexpo_clear = .
replace hbvexpo_clear = 1 if lower(anti_hbcii_interpretation)=="reactive" ///
                        & lower(hbsag_interpretation)=="nonreactive"
replace hbvexpo_clear = 0 if hbvexpo_clear==.
label values hbvexpo_clear YESNO
label var     hbvexpo_clear "Exposure & clearance (anti-HBc+ / HBsAg-)"

capture drop susceptible
gen byte susceptible = .
replace susceptible = 1 if lower(hbsag_interpretation)=="nonreactive" ///
                      & lower(anti_hbcii_interpretation)=="nonreactive" ///
                      & lower(hbsab_interpretation)=="nonreactive"
replace susceptible = 0 if susceptible==.
label values susceptible YESNO
label var     susceptible "Susceptible (HBsAg-/anti-HBc-/anti-HBs-)"

* Birth cohort (HBV vaccine epoch)
gen byte hbvdob = .
replace hbvdob = 0 if year(dateofbirth)<1995
replace hbvdob = 1 if inrange(year(dateofbirth),1995,1999)
replace hbvdob = 2 if inrange(year(dateofbirth),2000,2005)
label values hbvdob HBD
label var     hbvdob "HBV vaccine age epoch"

* HIV status & ART
gen byte hivstat = .
replace hivstat = 1 if hivelisa==1   // Positive
replace hivstat = 0 if hivelisa==2   // Negative
label values hivstat HIV
label var     hivstat "HIV status (ELISA)"

gen byte onart = .
replace onart = 0 if hivstat==1 & inrange(hivcascade,0,2)
replace onart = 1 if hivstat==1 & hivcascade>2 & hivcascade<.
label values onart YESNO
label var     onart "On ART (cascade)"

* HIV VL suppression (PLWH only; optional)
gen byte vl_suppressed = .
replace vl_suppressed = 1 if hivstat==1 & vl<50
replace vl_suppressed = 0 if hivstat==1 & vl>=50 & vl<.
label values vl_suppressed SUPP
label var     vl_suppressed "Viral load <50 copies/mL (PLWH only)"

* Hypertension & Diabetes (binary)
gen byte hypertension = .
replace hypertension = 1 if htncascade>0 & htncascade<.
replace hypertension = 0 if htncascade==0
label values hypertension YESNO
label var     hypertension "Hypertension"

gen byte diabetic = .
replace diabetic = 1 if dmdiag==1
replace diabetic = 0 if dmdiag==2
label values diabetic YESNO
label var     diabetic "Diabetes"

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

* S/CO threshold grid (strict)
foreach k in 1 2 5 10 100 {
    gen byte hbsag_gt`k' = (hbv_num>`k')  if hbv_num<.
    gen byte hbsag_le`k' = (hbv_num<=`k') if hbv_num<.
    label values hbsag_gt`k' YESNO
    label values hbsag_le`k' YESNO
    label var hbsag_gt`k' "HBsAg S/CO > `k' (strict)"
    label var hbsag_le`k' "HBsAg S/CO <= `k'"
}

* PRIMARY confirmed infection flag (S/CO >10)
gen byte hbv_active = (hbv_num>10) if hbv_num<.
label values hbv_active YESNO
label var     hbv_active "HBsAg S/CO >10 (confirmed infection; PRIMARY)"

compress
save "$DERIVED/vk_master.dta", replace

* Optional unweighted check (PI counts)
quietly count if hbsag_gt1==1
di as txt "Check S/CO>1 count: " as res r(N)
quietly count if hbsag_gt2==1
di as txt "Check S/CO>2 count: " as res r(N)
quietly count if hbsag_gt5==1
di as txt "Check S/CO>5 count: " as res r(N)
quietly count if hbsag_gt10==1
di as txt "Check S/CO>10 count: " as res r(N)
quietly count if hbsag_gt100==1
di as txt "Check S/CO>100 count: " as res r(N)

*-------------------------------*
* 2. Case–control lab merges -> hbv_cc_merged.dta
*-------------------------------*

* Cases (DNA/HBeAg/ALT panel)
import excel using "$DATA/EvoLVE AHRI_PLA_for HBV_VL-HBeAg-Alt_RESULTS 18 Oct 2024 1 Mod.xlsx", ///
    sheet("Sheet1") firstrow clear
rename _all, lower
tempfile hbv_cases
save `hbv_cases'

* Controls (ALT panel)
import excel using "$DATA/AHRI_SER_for Alt_RESULTS 2025-01-25 Mod.xlsx", ///
    sheet("Sheet1") firstrow clear
rename _all, lower
tempfile hbv_ctrl
save `hbv_ctrl'

* Merge each with VK master
use `hbv_cases', clear
merge 1:1 iintid using "$DERIVED/vk_master.dta", keep(match) nogen
gen byte case = 1
tempfile cases_m
save `cases_m'

use `hbv_ctrl', clear
merge 1:1 iintid using "$DERIVED/vk_master.dta", keep(match) nogen
gen byte case = 0
tempfile ctrls_m
save `ctrls_m'

* Append and label
use `ctrls_m', clear
append using `cases_m'
gen byte status = case
label values status CS
label var status "Case or control status"
compress
save "$DERIVED/hbv_cc_merged.dta", replace

*-------------------------------*
* 3. Harmonise DNA / ALT / HBeAg
*   and build canonical analysis file
*-------------------------------*

use "$DERIVED/hbv_cc_merged.dta", clear

* HBV DNA (IU/mL): bound/impute
gen str20 dna_vl_str = hbv_viral_load_iu_ml
replace dna_vl_str = "9.9"         if hbv_viral_load_iu_ml=="<10"
replace dna_vl_str = "1000000001"  if hbv_viral_load_iu_ml==">1000000000"
replace dna_vl_str = "0"           if hbv_viral_load_iu_ml=="Target Not Detected"
destring dna_vl_str, gen(dna_vl) force
label var dna_vl "HBV DNA IU/mL (clean; <10→9.9, TND→0, >1e9→1,000,000,001)"

gen double logvl = .
replace logvl = ln(dna_vl) if dna_vl>0
label var logvl "ln(HBV DNA IU/mL), 0 if TND handled pre-log"

* ALT (U/L): numeric
replace alt_u_l = subinstr(alt_u_l,"<7,00","6.99",.)
replace alt_u_l = subinstr(alt_u_l,"<7","6.99",.)
destring alt_u_l, replace force
label var alt_u_l "ALT (U/L)"

* ALT categories
gen byte ALT_category = .
replace ALT_category = 1 if (sex==1 & alt_u_l>30) | (sex==2 & alt_u_l>19)
replace ALT_category = 0 if missing(ALT_category) & alt_u_l<.
label values ALT_category ALTCAT
label var     ALT_category "ALT category (WHO ULN)"

gen byte nhls_altcat = .
replace nhls_altcat = 1 if (sex==1 & alt_u_l>40) | (sex==2 & alt_u_l>35)
replace nhls_altcat = 0 if missing(nhls_altcat) & alt_u_l<.
label values nhls_altcat ALTNHLS
label var     nhls_altcat "ALT category (NHLS ULN)"

* HBeAg (note: source var is 'hbeag_intepretation' in the sheet)
gen byte HBeAg_status = .
replace HBeAg_status = 1 if lower(hbeag_intepretation)=="reactive"
replace HBeAg_status = 0 if lower(hbeag_intepretation)=="nonreactive"
label define HBE 0 "HBeAg Negative" 1 "HBeAg Positive", replace
label values HBeAg_status HBE
label var     HBeAg_status "HBeAg status"

* Re-attach S/CO flags from master
preserve
use "$DERIVED/vk_master.dta", clear
keep iintid hbv_num hbv_active hbsag_gt1 hbsag_gt2 hbsag_gt5 hbsag_gt10 hbsag_gt100 ///
               hbsag_le1 hbsag_le2 hbsag_le5 hbsag_le10 hbsag_le100
tempfile flags
save `flags'
restore
merge 1:1 iintid using `flags', nogen keep(master match)

*-------------------------------*
* 3.5 Treatment eligibility flags (WHO-like)
*   PRIMARY: S/CO>10; Sensitivity: >5 and >2
*   Comorbidity = HIV not on ART OR diabetes OR hypertension 
*-------------------------------*

* Helper flags
gen byte dna_detect   = (dna_vl>=10) if dna_vl<.
gen byte alt_uln_elev = (ALT_category==1) if ALT_category<.
gen byte hiv_not_onart = (hivstat==1 & onart==0)
gen byte dm_yes = (diabetic==1)

* PRIMARY eligibility among S/CO>10
gen byte treateligible_primary = .
replace treateligible_primary = 0 if hbv_active==1
replace treateligible_primary = 1 if hbv_active==1 & ///
   ( (dna_vl>2000 & alt_uln_elev==1) | hiv_not_onart==1 | dm_yes==1 | hypertension==1 )
label values treateligible_primary YESNO
label var     treateligible_primary "Eligible (WHO-like) among S/CO>10 [PRIMARY]"

* Sensitivity: among S/CO>5
gen byte treateligible_gt5 = .
replace treateligible_gt5 = 0 if hbsag_gt5==1
replace treateligible_gt5 = 1 if hbsag_gt5==1 & ///
   ( (dna_vl>2000 & alt_uln_elev==1) | hiv_not_onart==1 | dm_yes==1 | hypertension==1 )
label values treateligible_gt5 YESNO
label var     treateligible_gt5 "Eligible (WHO-like) among S/CO>5 [Sensitivity]"

* Sensitivity: among S/CO>2
gen byte treateligible_gt2 = .
replace treateligible_gt2 = 0 if hbsag_gt2==1
replace treateligible_gt2 = 1 if hbsag_gt2==1 & ///
   ( (dna_vl>2000 & alt_uln_elev==1) | hiv_not_onart==1 | dm_yes==1 | hypertension==1 )
label values treateligible_gt2 YESNO
label var     treateligible_gt2 "Eligible (WHO-like) among S/CO>2 [Sensitivity]"

* HBV DNA category (descriptive)
gen byte hbvlcat = .
replace hbvlcat = 1 if dna_vl<10
replace hbvlcat = 2 if inrange(dna_vl,10,2000)
replace hbvlcat = 3 if inrange(dna_vl,2001,20000)
replace hbvlcat = 4 if dna_vl>20000
label values hbvlcat VLBIN
label var     hbvlcat "HBV DNA category (IU/mL)"

* Survey design
svyset [pweight=totwt]
keep if case==1 

compress
save "$DERIVED/analysis_eligibility.dta", replace

* QC prints 
di as res "Built: $DERIVED/analysis_eligibility.dta"
summ hbv_num if hbv_num<.
tab hbv_active if hbv_num<., m
svy: proportion treateligible_primary if hbv_active==1


*---------------------------------------------------------------*
* Load analysis file
*---------------------------------------------------------------*
use "$DERIVED/analysis_eligibility.dta", clear
svyset [pweight=totwt]

* Safety checks
foreach v in hbsag_gt1 hbsag_gt10 hivstat onart diabetic alt_u_l dna_vl sex {
    capture confirm variable `v'
    if _rc {
        di as err "Required variable `v' not found. Re-run builder."
        exit 198
    }
}

*---------------------------------------------------------------*
* Common flags
*---------------------------------------------------------------*
* WHO ALT ULN: >30 (men) / >19 (women)
cap drop alt_uln_elev
gen byte alt_uln_elev = ( (sex==1 & alt_u_l>30) | (sex==2 & alt_u_l>19) ) if alt_u_l<.

* DNA available (any numeric, incl TND=0)
cap drop dna_available
gen byte dna_available = (dna_vl < .)

* PLWH on ART (removed from evaluation as "already on treatment")
cap drop on_art_now
gen byte on_art_now = (hivstat==1 & onart==1)

* Comorbidities used for fallback eligibility
cap drop hiv_not_onart dm_yes
gen byte hiv_not_onart = (hivstat==1 & (onart==0 | onart==.))
gen byte dm_yes        = (diabetic==1)

*---------------------------------------------------------------*
* Helper that produces reason-coded counts for one threshold
*   FLOW("PRIMARY")  : pair first (ALT↑ & DNA>2000), else comorbidity
*   FLOW("ALT_ONLY") : ALT↑ first, else comorbidity
*---------------------------------------------------------------*
cap program drop __ppt_counts2
program define __ppt_counts2, rclass
    syntax, THRESHVAR(name) THRESHLAB(string) FLOW(string)

    tempvar pos notrt evalflag elig_pair not_pair elig_comorb final_inelig elig_alt not_alt
    gen byte `pos'   = (`threshvar'==1)
    gen byte `notrt' = (`pos'==1 & on_art_now==0)

    quietly count if `pos'==1
    return scalar n_pos    = r(N)

    quietly count if `pos'==1 & on_art_now==1
    return scalar n_onart  = r(N)

    quietly count if `notrt'==1
    return scalar n_notrt  = r(N)

    if "`flow'"=="PRIMARY" {
        * Evaluate only those with DNA available
        gen byte `evalflag' = (`notrt'==1 & dna_available==1)

        * Pair rule: ALT↑ & DNA>2000
        gen byte `elig_pair' = (`evalflag'==1 & alt_uln_elev==1 & dna_vl>2000)
        quietly count if `elig_pair'==1
        return scalar eligible_pair = r(N)

        * Those evaluated but not meeting the pair
        gen byte `not_pair' = (`evalflag'==1 & `elig_pair'==0)
        quietly count if `not_pair'==1
        return scalar not_pair = r(N)

        * Among those not_pair, eligible by comorbidity (HIV not on ART or DM)
        gen byte `elig_comorb' = (`not_pair'==1 & (hiv_not_onart==1 | dm_yes==1))
        quietly count if `elig_comorb'==1
        return scalar eligible_comorb = r(N)

        * Final ineligible = evaluated but neither pair nor comorbidity
        gen byte `final_inelig' = (`not_pair'==1 & `elig_comorb'==0)
        quietly count if `final_inelig'==1
        return scalar final_ineligible = r(N)

        * Evaluated counts + not evaluated (no DNA)
        quietly count if `evalflag'==1
        return scalar n_evaluated  = r(N)

        quietly count if `notrt'==1 & dna_available==0
        return scalar not_evaluated = r(N)

        return local flow "PRIMARY"
    }
    else if "`flow'"=="ALT_ONLY" {
        * Everyone not on treatment is evaluated
        gen byte `evalflag' = (`notrt'==1)

        * ALT-only rule
        gen byte `elig_alt' = (`evalflag'==1 & alt_uln_elev==1)
        quietly count if `elig_alt'==1
        return scalar eligible_alt = r(N)

        * Not ALT-eligible
        gen byte `not_alt' = (`evalflag'==1 & `elig_alt'==0)
        quietly count if `not_alt'==1
        return scalar not_alt = r(N)

        * Among not_alt, eligible by comorbidity (HIV not on ART or DM)
        gen byte `elig_comorb' = (`not_alt'==1 & (hiv_not_onart==1 | dm_yes==1))
        quietly count if `elig_comorb'==1
        return scalar eligible_comorb = r(N)

        * Final ineligible = not_alt but no comorbidity
        gen byte `final_inelig' = (`not_alt'==1 & `elig_comorb'==0)
        quietly count if `final_inelig'==1
        return scalar final_ineligible = r(N)

        * Evaluated (all) and not evaluated (none)
        quietly count if `evalflag'==1
        return scalar n_evaluated  = r(N)
        return scalar not_evaluated = 0

        return local flow "ALT_ONLY"
    }
    else {
        di as err "FLOW must be PRIMARY or ALT_ONLY"
        exit 198
    }

    return local threshold "`threshlab'"
end

*---------------------------------------------------------------*
* Run for >1 and >10; export reason-coded columns for PPT
*---------------------------------------------------------------*
tempname H
postutil clear
tempfile T
postfile `H' str10 flow str5 threshold ///
    double n_pos n_onart n_notrt n_evaluated not_evaluated ///
           eligible_pair not_pair eligible_alt not_alt ///
           eligible_comorb final_ineligible ///
    using "`T'", replace

* PRIMARY flow (pair first)
quietly __ppt_counts2, threshvar(hbsag_gt1)  threshlab(">1")  flow(PRIMARY)
post `H' ("PRIMARY") (">1")  (r(n_pos)) (r(n_onart)) (r(n_notrt)) (r(n_evaluated)) (r(not_evaluated)) ///
                      (r(eligible_pair)) (r(not_pair)) (.) (.) (r(eligible_comorb)) (r(final_ineligible))

quietly __ppt_counts2, threshvar(hbsag_gt10) threshlab(">10") flow(PRIMARY)
post `H' ("PRIMARY") (">10") (r(n_pos)) (r(n_onart)) (r(n_notrt)) (r(n_evaluated)) (r(not_evaluated)) ///
                      (r(eligible_pair)) (r(not_pair)) (.) (.) (r(eligible_comorb)) (r(final_ineligible))

* ALT-ONLY flow
quietly __ppt_counts2, threshvar(hbsag_gt1)  threshlab(">1")  flow(ALT_ONLY)
post `H' ("ALT_ONLY") (">1")  (r(n_pos)) (r(n_onart)) (r(n_notrt)) (r(n_evaluated)) (r(not_evaluated)) ///
                      (.) (.) (r(eligible_alt)) (r(not_alt)) (r(eligible_comorb)) (r(final_ineligible))

quietly __ppt_counts2, threshvar(hbsag_gt10) threshlab(">10") flow(ALT_ONLY)
post `H' ("ALT_ONLY") (">10") (r(n_pos)) (r(n_onart)) (r(n_notrt)) (r(n_evaluated)) (r(not_evaluated)) ///
                      (.) (.) (r(eligible_alt)) (r(not_alt)) (r(eligible_comorb)) (r(final_ineligible))

postclose `H'

use "`T'", clear
order flow threshold n_pos n_onart n_notrt n_evaluated not_evaluated ///
      eligible_pair not_pair eligible_alt not_alt eligible_comorb final_ineligible
format n_pos-final_ineligible %9.0gc

label var flow              "Decision flow"
label var threshold         "HBsAg threshold"
label var n_pos             "HBsAg positives (strict >k)"
label var n_onart           "Already on treatment (PLWH on ART)"
label var n_notrt           "Not on treatment"
label var n_evaluated       "Evaluated"
label var not_evaluated     "Not evaluated (no DNA)"
label var eligible_pair     "Eligible: ALT↑ & DNA>2000"
label var not_pair          "Not eligible by pair"
label var eligible_alt      "Eligible: ALT↑ (ALT-only flow)"
label var not_alt           "Not eligible by ALT (ALT-only flow)"
label var eligible_comorb   "Eligible by comorbidity (HIV+ no ART or diabetes)"
label var final_ineligible  "Final NOT eligible"

list, sepby(flow) noobs abbreviate(20)

export excel using "$OUT/hbv_treatment_eligibility_counts.xlsx", ///
    firstrow(variables) replace
di as res "Exported: $OUT/hbv_treatment_eligibility_counts.xlsx"
