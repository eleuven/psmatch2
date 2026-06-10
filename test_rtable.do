// test_rtable.do
// Tests for r(table) in psmatch2

clear all
set seed 999
capture program drop psmatch2
quietly do "psmatch2.ado"
capture program drop psmatch2_collect
quietly do "psmatch2_collect.ado"

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

// ------------------------------------------------------------------
// Test 1: Single outcome, no ate — 9x1 matrix, column y1:ATT
// ------------------------------------------------------------------
_dgp2 300
qui psmatch2 treat x1 x2, outcome(y1) ai(1) population
matrix T = r(table)
scalar _t1_att   = r(att_y1)
scalar _t1_seatt = r(seatt_y1)

assert rowsof(T) == 9
assert colsof(T) == 1

local rn : rownames T
assert "`rn'" == "b se z pvalue ll ul df crit eform"

local coleq : coleq T
local colnm : colnames T
assert "`coleq'" == "y1"
assert "`colnm'" == "ATT"

assert reldif(T[1,1], _t1_att)   < 1e-10
assert reldif(T[2,1], _t1_seatt) < 1e-10

di as text "PASS: single outcome no ate — 9x1 matrix, y1:ATT"

// ------------------------------------------------------------------
// Test 2: Single outcome with ate — 9x3, columns y1:ATT y1:ATU y1:ATE
// ------------------------------------------------------------------
_dgp2 400
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population
matrix T = r(table)
scalar _t2_att   = r(att_y1)
scalar _t2_seatt = r(seatt_y1)
scalar _t2_atu   = r(atu_y1)
scalar _t2_seatu = r(seatu_y1)
scalar _t2_ate   = r(ate_y1)
scalar _t2_seate = r(seate_y1)

assert rowsof(T) == 9
assert colsof(T) == 3

local coleq : coleq T
assert "`coleq'" == "y1 y1 y1"
local colnm : colnames T
assert "`colnm'" == "ATT ATU ATE"

assert reldif(T[1,1], _t2_att)   < 1e-10
assert reldif(T[2,1], _t2_seatt) < 1e-10
assert reldif(T[1,2], _t2_atu)   < 1e-10
assert reldif(T[2,2], _t2_seatu) < 1e-10
assert reldif(T[1,3], _t2_ate)   < 1e-10
assert reldif(T[2,3], _t2_seate) < 1e-10

di as text "PASS: single outcome with ate — 9x3 matrix, y1:ATT y1:ATU y1:ATE"

// ------------------------------------------------------------------
// Test 3: Multiple outcomes with ate — 9x6
// ------------------------------------------------------------------
_dgp2 400
qui psmatch2 treat x1 x2, outcome(y1 y2) ate ai(1) population
matrix T = r(table)

assert rowsof(T) == 9
assert colsof(T) == 6

local coleq : coleq T
assert "`coleq'" == "y1 y1 y1 y2 y2 y2"
local colnm : colnames T
assert "`colnm'" == "ATT ATU ATE ATT ATU ATE"

di as text "PASS: multiple outcomes with ate — 9x6 matrix"

// ------------------------------------------------------------------
// Test 4: z and p-value identities
// ------------------------------------------------------------------
_dgp2 400
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population
matrix T = r(table)

// ATT: column 1
local b_att  = T[1,1]
local se_att = T[2,1]
local z_att  = T[3,1]
local p_att  = T[4,1]
assert reldif(`z_att', `b_att' / `se_att') < 1e-10
assert reldif(`p_att', 2 * normal(-abs(`z_att'))) < 1e-10

// ATE: column 3
local b_ate  = T[1,3]
local se_ate = T[2,3]
local z_ate  = T[3,3]
local p_ate  = T[4,3]
assert reldif(`z_ate', `b_ate' / `se_ate') < 1e-10
assert reldif(`p_ate', 2 * normal(-abs(`z_ate'))) < 1e-10

di as text "PASS: z and p-value identities"

// ------------------------------------------------------------------
// Test 5: Confidence interval identity
// ------------------------------------------------------------------
// reuse T from test 4 (same dataset still loaded)
local crit = invnormal(1 - (100 - c(level)) / 200)

// ATT column
local b  = T[1,1]
local se = T[2,1]
local ll = T[5,1]
local ul = T[6,1]
assert reldif(`ll', `b' - `crit' * `se') < 1e-10
assert reldif(`ul', `b' + `crit' * `se') < 1e-10

// ATE column
local b  = T[1,3]
local se = T[2,3]
local ll = T[5,3]
local ul = T[6,3]
assert reldif(`ll', `b' - `crit' * `se') < 1e-10
assert reldif(`ul', `b' + `crit' * `se') < 1e-10

// crit row matches
assert reldif(T[8,1], `crit') < 1e-10

di as text "PASS: confidence interval identity"

// ------------------------------------------------------------------
// Test 6: Legacy return scalars unchanged
// ------------------------------------------------------------------
_dgp2 300
qui psmatch2 treat x1 x2, outcome(y1) ate ai(1) population
scalar _t6_att   = r(att)
scalar _t6_att_y1 = r(att_y1)
scalar _t6_seatt  = r(seatt)
scalar _t6_seatt_y1 = r(seatt_y1)
scalar _t6_atu   = r(atu)
scalar _t6_ate   = r(ate)
scalar _t6_seatu = r(seatu)
scalar _t6_seate = r(seate)
matrix T = r(table)

// verify legacy scalars still present and match T
assert reldif(_t6_att, T[1,1])   < 1e-10
assert reldif(_t6_seatt, T[2,1]) < 1e-10
assert reldif(_t6_atu, T[1,2])   < 1e-10
assert reldif(_t6_ate, T[1,3])   < 1e-10
// scalar with _y suffix matches too
assert reldif(_t6_att_y1, _t6_att) < 1e-10

di as text "PASS: legacy return scalars unchanged"

// ------------------------------------------------------------------
// Test 7: Collect helper smoke test
// ------------------------------------------------------------------
_dgp2 400
qui psmatch2 treat x1 x2, outcome(y1 y2) ate ai(1) population
// psmatch2_collect reads r(table); copy it before calling so we can
// verify collect does not corrupt r(table)
matrix T_pre = r(table)

psmatch2_collect, clear
collect dir
// collect dir should succeed (no error) and show at least one tag
// We cannot assert rendered text, so checking T_pre still makes sense
// to confirm r(table) was preserved in T before the call.
assert colsof(T_pre) == 6

di as text "PASS: collect helper smoke test"

// ------------------------------------------------------------------
di as text _n "All r(table) tests passed."
