// test_ai2016.do
// Tests for AI(2016) estimated-propensity-score SE correction in psmatch2

clear all
set seed 42
qui do "/Users/edwinl/Documents/GitHub/psmatch2/psmatch2.ado"

// Two continuous covariates, probit treatment model, two outcomes
capture program drop _dgp2
program define _dgp2
    args n
    clear
    set obs `n'
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen double idx = 0.4*x1 + 0.6*x2
    gen treat = (rnormal() < idx)
    gen y1 = x1 + x2 + treat + rnormal()
    gen y2 = 2*x1 - x2 + 0.5*treat + rnormal()
end


// -------------------------------------------------------------------------
// Test 1: Backward compatibility — ai(0) results identical across two runs
// -------------------------------------------------------------------------

di as text _n "=== Test 1: Backward compatibility (ai(0)) ==="
set seed 101
_dgp2 300

qui psmatch2 treat x1 x2, outcome(y1) ate ai(0) population
scalar _t1_att   = r(att)
scalar _t1_seatt = r(seatt)
scalar _t1_ate   = r(ate)
scalar _t1_seate = r(seate)
local  _t1_qa    = r(qA_y1)

qui psmatch2 treat x1 x2, outcome(y1) ate ai(0) population
assert reldif(r(att),   _t1_att)   < 1e-10
assert reldif(r(seatt), _t1_seatt) < 1e-10
assert reldif(r(ate),   _t1_ate)   < 1e-10
assert reldif(r(seate), _t1_seate) < 1e-10
assert missing(`_t1_qa')
assert missing(r(qA_y1))

di as text "  PASS: ai(0) backward compatibility"


// -------------------------------------------------------------------------
// Test 2: Matching output preserved when correction fires
// Internal-PS eligible run vs external-PS reference (same probit fit)
// -------------------------------------------------------------------------

di as text _n "=== Test 2: Matching output preserved ==="
set seed 202
_dgp2 400

// eligible run (internal PS estimation)
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population
scalar _t2_att = r(att)
scalar _t2_atu = r(atu)
scalar _t2_ate = r(ate)
tempvar n1_ref ps_ref wt_ref
gen long   `n1_ref' = _n1
gen double `ps_ref' = _pscore
gen double `wt_ref' = _weight

// reference run: same probit, external pscore makes correction ineligible
qui probit treat x1 x2
qui predict double _ps_ext2, pr
qui psmatch2 treat, pscore(_ps_ext2) outcome(y1) ate ai(1) population
scalar _t2_att_ref = r(att)
scalar _t2_atu_ref = r(atu)
scalar _t2_ate_ref = r(ate)
drop _ps_ext2

// matching output identical
gen byte _t2_ok_ps = abs(_pscore - `ps_ref') < 1e-10
qui sum _t2_ok_ps
assert r(min) == 1
drop _t2_ok_ps

gen byte _t2_ok_n1 = (_n1 == `n1_ref') | (missing(_n1) & missing(`n1_ref'))
qui sum _t2_ok_n1
assert r(min) == 1
drop _t2_ok_n1

gen byte _t2_ok_wt = (abs(_weight - `wt_ref') < 1e-10) | (missing(_weight) & missing(`wt_ref'))
qui sum _t2_ok_wt
assert r(min) == 1
drop _t2_ok_wt

// point estimates identical
assert reldif(_t2_att_ref, _t2_att) < 1e-10
assert reldif(_t2_atu_ref, _t2_atu) < 1e-10
assert reldif(_t2_ate_ref, _t2_ate) < 1e-10

di as text "  PASS: matching output unchanged when correction fires"


// -------------------------------------------------------------------------
// Test 3: Ineligible — pscore() supplied externally
// -------------------------------------------------------------------------

di as text _n "=== Test 3: Ineligible (pscore() supplied) ==="
set seed 303
_dgp2 300

qui probit treat x1 x2
qui predict double _ps_ext3, pr

qui psmatch2 treat, pscore(_ps_ext3) outcome(y1) ate ai(1) population
assert missing(r(qA_y1))
scalar _t3_seate = r(seate)

qui psmatch2 treat, pscore(_ps_ext3) outcome(y1) ate ai(1) population
assert reldif(r(seate), _t3_seate) < 1e-10
drop _ps_ext3

di as text "  PASS: correction not applied when pscore() supplied"


// -------------------------------------------------------------------------
// Test 4: Eligible probit — qA >= 0, seate corrected, point estimates unchanged
// seate_ai_fixed_y1 == seate from external-pscore reference run
// -------------------------------------------------------------------------

di as text _n "=== Test 4: Eligible probit ==="
set seed 404
_dgp2 500

// reference: external pscore run, correction ineligible
qui probit treat x1 x2
qui predict double _ps_ext4, pr
qui psmatch2 treat, pscore(_ps_ext4) outcome(y1) ate ai(1) population
scalar _t4_seate_ref = r(seate)
scalar _t4_att_ref   = r(att)
scalar _t4_atu_ref   = r(atu)
scalar _t4_ate_ref   = r(ate)
drop _ps_ext4

// eligible run
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population
scalar _t4_qA       = r(qA_y1)
scalar _t4_seate    = r(seate)
scalar _t4_seate_fx = r(seate_ai_fixed_y1)

// qA is a quadratic form: non-negative
assert _t4_qA >= -1e-10
// seate_ai_fixed matches pre-correction reference
assert reldif(_t4_seate_fx, _t4_seate_ref) < 1e-10
// corrected SE weakly below uncorrected
assert _t4_seate <= _t4_seate_fx + 1e-8
// point estimates unchanged
assert reldif(r(att), _t4_att_ref) < 1e-10
assert reldif(r(atu), _t4_atu_ref) < 1e-10
assert reldif(r(ate), _t4_ate_ref) < 1e-10

di as text "  PASS: probit eligible — qA>=0, seate corrected, point estimates unchanged"


// -------------------------------------------------------------------------
// Test 5: SE identity — seate^2 = seate_ai_fixed^2 - qA
// -------------------------------------------------------------------------

di as text _n "=== Test 5: SE identity ==="
// data still in memory from Test 4
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population

if (r(seate) < .) {
    assert abs(r(seate)^2 - (r(seate_ai_fixed_y1)^2 - r(qA_y1))) < 1e-10
}
if (r(seatt) < .) {
    assert abs(r(seatt)^2 - (r(seatt_ai_fixed_y1)^2 - r(qTminus_y1) + r(qTplus_y1))) < 1e-10
}
if (r(seatu) < .) {
    assert abs(r(seatu)^2 - (r(seatu_ai_fixed_y1)^2 - r(qUminus_y1) + r(qUplus_y1))) < 1e-10
}

di as text "  PASS: SE identity holds"


// -------------------------------------------------------------------------
// Test 6: Eligible logit
// -------------------------------------------------------------------------

di as text _n "=== Test 6: Eligible logit ==="
set seed 606
_dgp2 400

qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population logit
assert r(qA_y1) >= -1e-10
assert r(qA_y1) < .
assert r(seate) <= r(seate_ai_fixed_y1) + 1e-8

di as text "  PASS: logit eligible — correction fires and qA>=0"


// -------------------------------------------------------------------------
// Test 7: Multiple outcomes — corrections are outcome-specific
// -------------------------------------------------------------------------

di as text _n "=== Test 7: Multiple outcomes ==="
set seed 707
_dgp2 400

qui psmatch2 treat x1 x2, outcome(y1 y2) ate ai(1) population
assert r(qA_y1) >= -1e-10
assert r(qA_y2) >= -1e-10
assert r(qA_y1) < .
assert r(qA_y2) < .
assert r(seate_y1) < .
assert r(seate_y2) < .
assert r(seate_ai_fixed_y1) < .
assert r(seate_ai_fixed_y2) < .

di as text "  PASS: multiple outcomes — corrections outcome-specific"


// -------------------------------------------------------------------------
// Test 8: Factor variables — correction skipped, psmatch2 runs without error
// -------------------------------------------------------------------------

di as text _n "=== Test 8: Factor variables ==="
set seed 808
_dgp2 300
gen x_cat = ceil(3 * runiform())

qui psmatch2 treat i.x_cat x2, outcome(y1) ate ai(1) population
assert missing(r(qA_y1))
assert r(ate) < .

di as text "  PASS: factor variables — psmatch2 runs, correction skipped"


// -------------------------------------------------------------------------
// Test 9: Other ineligible conditions — no error, correction skipped
// -------------------------------------------------------------------------

di as text _n "=== Test 9: Ineligible conditions ==="
set seed 909
_dgp2 300

// 9a: positive caliper
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population caliper(0.5)
assert missing(r(qA_y1))

// 9b: ties
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population ties
assert missing(r(qA_y1))

// 9c: noreplacement
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population noreplacement
assert missing(r(qA_y1))

// 9d: altvariance
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population altvariance
assert missing(r(qA_y1))

// 9e: no population
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1)
assert missing(r(qA_y1))

// 9f: no ate
qui psmatch2 treat x1 x2, outcome(y1) ai(1) population
assert missing(r(qA_y1))

di as text "  PASS: ineligible conditions — no error, correction skipped"


di as text _n "All tests passed."
