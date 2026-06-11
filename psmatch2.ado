*! version 4.0.13 30may2026 E. Leuven, B. Sianesi
program define psmatch2, sortpreserve
	version 11.0
	#delimit ;
	syntax varlist(min=1 fv) [if] [in] [,
	OUTcome(varlist)
	Pscore(varname)
	Neighbor(integer 1)
	TIES
	RADIUS
	CALiper(real 0)
	MAHALanobis(varlist fv)
	KERNEL
	LLR
	Kerneltype(string)
	BWidth(string)
	COMmon
	AI(integer 0)
	POPulation
	ALTVariance
	TRIM(real 100)
	ODDS
	LOGIT
	INDEX
	QUIetly
	NOREPLacement
	DESCending
	WARNings
	ATE
	W(string)
	SPLINE
	NKnots(integer 0)
	];
	#delimit cr

	// record sort order
	tempvar order
	g long `order' = _n

	// clean up data
	foreach v in _treated _support _weight _pscore _id _nn _pdif {
		cap drop `v'
	}
	forv n=1 / `neighbor' {
		cap drop _n`n'
	}
	global OUTVAR `outcome'
	if ("`outcome'"!="") {
		foreach v of varlist `outcome' {
			cap drop _`v'
			local moutvar `moutvar' _`v'
		}
	}

	// determine subset we work on
	marksample touse
	capture markout `touse' `outcome' `control' `mahalanobis'

	// separate treatment indicator from varlist
	tokenize `varlist'
	local treat `1'
	macro shift
	local varlist "`*'"
	local k : word count `varlist'

	// determine matching metric
	if ("`mahalanobis'"!="") local metric = "mahalanobis"
	else local metric = "pscore"

	if (`k'==0 & "`pscore'"=="" & "`metric'"=="pscore") {
		di as error "You should either specify a " as input "varlist" as error " or " as input "propensity score"
		exit 198
	}
	if (`k'>0 & "`pscore'"!="") {
		di as error "You cannot specify both a " as input "varlist" as error " AND a " as input "propensity score"
		exit 198
	}
	
	if ("`ties'"!="" & "`metric'"=="mahalanobis") {
		di as text "Ties is not implemented for Mahalanobis matching. This option will be ignored."
	}

	// determine matching method
	local method "neighbor"
	if ("`kernel'"!="") local method "`kernel'"
	if ("`llr'"!="") local method "`llr'"

	if ("`llr'"!="" & "`kernel'"!="") {
		di as error "You cannot do kernel and llr matching at the same time"
		exit 198
	}

	if ("`noreplacement'"!="" & ("`method'"!="neighbor" | "`metric'"!="pscore" | `neighbor'>1 )) {
		di as error "Matching without replacement is only implemented with 1-to-1 propensity score matching"
		exit 198
	}
	if ("`descending'"!="" & "`noreplacement'"=="") {
		di as error "Option " as input "descending" as error " makes only sense when matching without replacement"
		exit 198
	}

	//
	if (`ai' > 0 & !("`method'" == "neighbor")) {
		di as error "Option -ai- only allowed when doing nearest neighbor matching on the covariates (Mahalanobis)"
		exit 198
	}

	// set caliper to missing (infinity) if not requested
	local caliper_orig = `caliper'
	if (`caliper'==0) local caliper = .

	// check kerneltype
	if ("`method'"=="kernel" & "`kerneltype'"=="") local kerneltype "epan"
	if ("`method'"=="llr" & "`kerneltype'"=="") local kerneltype "epan"

	if !("`kerneltype'"=="" | "`kerneltype'"=="normal" | "`kerneltype'"=="epan" | "`kerneltype'"=="biweight" | "`kerneltype'"=="uniform" | "`kerneltype'"=="tricube") {
		di as error "Kerneltype `kerneltype' not recognized"
		exit 198
	}

	// radius matching is like kernel matching with uniform kernel
	if ("`radius'"!="") {
		local method "kernel"
		local kerneltype "uniform"
		local bwidth = `caliper'
	}
	if ("`bwidth'"=="") local bwidth "0.06"

	// AI(2016): internal pscore NN, population AI
	local do_pscorr = (`ai' > 0                                      ///
		& "`population'" != ""                                       ///
		& "`ate'" != ""                                              ///
		& "`method'" == "neighbor"                                   ///
		& "`metric'" == "pscore"                                     ///
		& "`varlist'" != ""                                          ///
		& "`pscore'" == ""                                           ///
		& "`index'`odds'`ties'`noreplacement'`altvariance'`common'" == "" ///
		& `caliper_orig' <= 0                                        ///
		& `trim' == 100                                              ///
		& "`kernel'`llr'`radius'`spline'`mahalanobis'" == "")

	local psfix_note = (`ai' > 0 & "`metric'" == "pscore" & !`do_pscorr')

	// estimate propensity score
	if ("`varlist'"!="") {
		if ("`logit'"=="") {
			local logit "probit"
		}
		`quietly' `logit' `treat' `varlist' if `touse', nolog
		qui replace `touse' = e(sample) // factor vars may drop obs
		tempvar pscore
		qui predict double `pscore', `index'
		qui g double _pscore = `pscore' if `touse'
		label var _pscore "psmatch2: Propensity Score"
		// first-stage objects for AI(2016)
		if (`do_pscorr') {
			tempname Vgamma
			matrix `Vgamma' = e(V)
			local gammacols : colnames `Vgamma'
			local psxvars ""
			foreach _gc of local gammacols {
				if "`_gc'" != "_cons" local psxvars `psxvars' `_gc'
			}
			foreach x of local psxvars {
				capture confirm variable `x'
				if _rc {
					local do_pscorr = 0
					local psfix_note = 1
				}
			}
			if (`do_pscorr') {
				tempvar xb_ps dP_ps
				qui predict double `xb_ps', xb
				if "`logit'" == "logit" {
					qui gen double `dP_ps' = _pscore * (1 - _pscore) if `touse'
				}
				else {
					qui gen double `dP_ps' = normalden(`xb_ps') if `touse'
				}
			}
		}
	}
	else if ("`metric'"=="pscore") {
		qui g double _pscore = `pscore' 
		label var _pscore "psmatch2: Propensity Score"
	}
	capture markout `touse' _pscore

	// match on log odds ratio if requested, only with logit?
	if ("`odds'"!="") qui replace _pscore = ln(_pscore/(1 - _pscore))

	// create treatment indicator variable
	qui g byte _treated = `treat' if `touse'
	label variable _treated "psmatch2: Treatment assignment"
	cap label drop _treated
	label define _treated 0 "Untreated" 1 "Treated"
	label value _treated _treated

	// common support if requested
	if (("`common'"!="" | `trim'<100) & ("`varlist'"=="" & "`pscore'"=="")) {
		di as error "With option 'common' a propensity score is needed. Provide one"
		di as error "with option 'pscore()' or estimate one. See the help file for more details."
		exit 198
	}
	qui g byte _support = 1 if `touse'
	label variable _support "psmatch2: Common support"
	cap label drop _support
	label define _support 0 "Off support" 1 "On support"
	label value _support _support
	if (("`common'"!="" | `trim'<100) & ("`varlist'"!="" | "`pscore'"!="")) {
		if !inrange(`trim',0,100) {
			di as error "Trim level out of range"
			exit 198
		}
		qui _Support_ `pscore', level(`trim') `ate'
	}

	// do nearest neighbor if llr with tricube
	if ("`method'"=="llr" & "`kerneltype'"=="epan" & "`metric'"=="pscore") {
		local method "neighbor"
		if ("`bwidth'"!="") {
			local bwidth bw(`bwidth')
		}
		global OUTVAR ""
		foreach v of varlist `outcome' {
			cap drop _s_`v'
			qui lpoly `v' _pscore if _treated==0 & _support==1, nograph deg(1) at(_pscore) gen(_s_`v') `bwidth'
			if ("`ate'"!="") {
				tempvar s`v'
				qui lpoly `v' _pscore if _treated==1 & _support==1, nograph deg(1) at(_pscore) gen(`s`v'') `bwidth'
				qui replace _s_`v' = `s`v'' if _treated==1 & _support==1
			}
			global OUTVAR $OUTVAR _s_`v'
			label var _s_`v' "psmatch2: smoothed outcome variable"
		}
	}

	// spline
	if ("`spline'"!="") {
		cap which spline
		if (_rc) {
			di as error "You need to install package -spline-, for option spline"
			exit  198
		}
		local method "neighbor"
		if ("`nknots'"=="0") {
			qui count if _treated==0 & _support==1
			local nknots = int(r(N)^0.25)
		}
		global OUTVAR ""
		foreach v of varlist `outcome' {
			cap drop _s_`v'
			qui spline `v' _pscore if _treated==0 & _support==1, gen(_s_`v') nknots(`nknots') nograph
			if ("`ate'"!="") {
				if ("`nknots'"=="0") {
					qui count if _treated==1 & _support==1
					local nknots = int(r(N)^0.25)
				}
				tempvar s`v'
				qui spline `v' _pscore if _treated==1 & _support==1, gen(`s`v'') nknots(`nknots') nograph
				qui replace _s_`v' = `s`v'' if _treated==1 & _support==1
			}
			global OUTVAR $OUTVAR _s_`v'
			label var _s_`v' "psmatch2: smoothed outcome variable using -spline-"
		}
	}

	// create vars we will need
	qui g double _weight = _treated if _support==1
	char _weight[Type] "aweight"
	if "`ate'"!="" {
		qui replace _weight = 0 if _treated==1 & _support==1
	}
	label var _weight "psmatch2: weight of matched controls"

	// outcome of matches
	if ("`outcome'"!="") {
		foreach v of varlist `outcome'	 {
			if ("`ate'"=="") {
				qui g double _`v' = 0 if _support==1 & _treated==1
			}
			else qui g double _`v' = 0 if _support==1
			label var _`v' "psmatch2: value of `v' of match(es)"
		}
	}

	// check for duplicate pscores
	if ("`warnings'"!="" & "`metric'"=="pscore") {
		sort _treated _pscore
		cap by _treated _pscore: assert _N==1 if _treated==0 & _support==1
		if (!_rc & "`ate'"!="") {
			cap by _treated _pscore: assert _N==1 if _treated==1 & _support==1
		}
		if (_rc & "`method'"=="neighbor") {
			di as res "There are observations with identical propensity score values."
			di as res "The sort order of the data could affect your results."
			di as res "Make sure that the sort order is random before calling psmatch2."
		}
	}

	// sort data on treatment status and pscore and create id
	tempvar msup
	qui g byte `msup' = - _support
	if ("`metric'"=="pscore" & "`method'"=="neighbor") {
		if ("`descending'"=="") {
			sort `msup' _treated _pscore `order'
		}
		else {
			tempvar mps
			qui g double `mps' = - _pscore
			sort `msup' _treated `mps' `order'
		}
	}
	else sort `msup' _treated `order'

	if ("`method'"=="neighbor") {
		g _id = _n
		label var _id "psmatch2: Identifier (ID)"
		qui compress _id
		local idtype : type _id
		forv n=1/`neighbor' {
			qui g `idtype' _n`n' = .
			label var _n`n' "psmatch2: ID of nearest neighbor nr. `n'"
		}
		qui g _nn = 0 if _support==1
		label var _nn "psmatch2: # matched neighbors"
	}
	
	if ("`mahalanobis'" != "") {
		_rmdcoll _treated `mahalanobis' if `touse', expand
		local mahalanobis `r(varlist)'
	}	
	
	// calculate within sample covariance matrix if necessary
	if ("`metric'"=="mahalanobis") {
		if ("`w'"=="") {
			tempname XX XX0 XX1 w
			qui mat accum `XX'  = `mahalanobis' if _treated<=1, dev noc //XX
			local mahalanobis : rowfullnames `XX' //XX
			qui mat accum `XX0' = `mahalanobis' if _treated==0, dev noc
			qui mat accum `XX1' = `mahalanobis' if _treated==1, dev noc
			qui count if _treated<=1 & _support==1
			mat `w' = syminv((`XX0' + `XX1')/(r(N) - 2))
		}
		local matchon `mahalanobis'
	}
	else local matchon `pscore'

	/*
Your -st_addvar()- specification actually creates a variable of type
float. It is just that the variable index is not what you expected. The
problem is this: if there are estimation results stored in -e()- that
contain the -e(sample)- macro, there will be a hidden variable in the
dataset that stores an indicator (of type byte) that indicates the
estimation sample. After an estimation command terminates, this hidden
variable is stored at the last position of the dataset. Now, what
happened prior to Stata 15 is that using -st_addvar()- could enable
users to create new variables after the -e(sample)- variable and that
caused a number of problems in different contexts. In Stata 15, we made
a change to this behavior such that whenever a new view is created, the
-e(sample)- variable is automatically being moved to the end of the
dataset. In your case, this new behavior has the unfortunate side effect
that the index you are creating for the new temporary variable is now
referring to the estimation sample indicator instead of your new
variable (hence -st_vartype()- returns -byte- instead of -float-).

I think the most straightforward workaround to this would be to do a
-ereturn clear- once you have grabbed the -e()- results you need prior
to performing the matching. In that case there would be no -e(sample)-
variable in the dataset and everything should work like before.
		*/
	ereturn clear

	if ("`method'"=="neighbor" ) {
		qui count if _treated<=1 & _support==1
		local N = r(N)
		qui count if _treated==0 & _support==1
		local N0 = r(N)

		unab n1 : _n1-_n`neighbor'
	
		if ("`ties'"!="") local ties 1
		else local ties 0

		if ("`noreplacement'"!="") local noreplace 1
		else local noreplace 0

		if ("`ate'"!="") {
			mata : match_`metric'(1, `N0', `=`N0'+1', `N', `neighbor', `caliper', `noreplace', `ties', "`w'", "`mahalanobis'", "`n1'", "$OUTVAR", "`moutvar'", "`ate'")
			qui replace _support = 0 if _n1>=. & _treated==0
		}
		mata : match_`metric'(`=`N0'+1', `N', 1, `N0', `neighbor', `caliper', `noreplace', `ties', "`w'", "`mahalanobis'", "`n1'", "$OUTVAR", "`moutvar'", "`ate'")
		qui replace _support = 0 if _n1>=. & _treated==1

		// difference pscore between treat obs and nearest match
		if ("`metric'"=="pscore") {
			qui g double _pdif = abs(_pscore - _pscore[_n1])
			label var _pdif "psmatch2: abs(pscore - pscore[nearest neighbor])"
		}
		if (`ai' > 0) {
			// outcome of matches, treated to treated, controls to controls
			if ("`outcome'"!="") {
				foreach v of varlist `outcome'	 {
					cap drop _self_`v'
					qui g double _self_`v' = 0 if _support==1
					label var _self_`v' "psmatch2: matched value of `v' (T-T & C-C)"
					local soutvar `soutvar' _self_`v'
				}
			}
			// same-arm matched X means for AI(2016)
			local self_outvar  "$OUTVAR"
			local self_moutvar "`soutvar'"
			local selfxvars ""
			if (`do_pscorr') {
				local xi = 0
				foreach x of local psxvars {
					local ++xi
					tempvar selfx`xi'
					qui gen double `selfx`xi'' = 0 if _support==1
					local selfxvars `selfxvars' `selfx`xi''
				}
				local self_outvar  "`self_outvar' `psxvars'"
				local self_moutvar "`self_moutvar' `selfxvars'"
			}
			// match controls to controls
			mata : match_`metric'(1, `N0', 1, `N0', `ai', `caliper', `noreplace', `ties', "`w'", "`mahalanobis'", "`n1'", "`self_outvar'", "`self_moutvar'", "`ate'")
			// match treated to treated
			mata : match_`metric'(`=`N0' + 1', `N', `=`N0' + 1', `N', `ai', `caliper', `noreplace', `ties', "`w'", "`mahalanobis'", "`n1'", "`self_outvar'", "`self_moutvar'", "`ate'")
		}
	}
	else { // llr and kernel
		qui _Match_`method' `matchon', out(`outcome') metric(`metric') kerneltype(`kerneltype') bw(`bwidth') w(`w') `ate'
	}

	// controls off support
	qui replace _weight = . if _weight==0 | _support==0

	// generate output
	local _psc_obj ""
	if (`do_pscorr') {
		local _psc_obj "dp(`dP_ps') xvars(`psxvars') selfxvars(`selfxvars') vgamma(`Vgamma')"
	}
	_mktab `outcome', `ate' `spline' `llr' k(`kerneltype') ai(`ai') n(`neighbor') ///
		`population' `altvariance' exog(`varlist') ///
		pscorr(`do_pscorr') psfixnote(`psfix_note') `_psc_obj'

	// get rid of evil global
	macro drop OUTVAR
end



// FORMAT OUTPUT TABLE
program define _mktab, rclass
syntax [varlist(default=none)] [, ate spline llr Kerneltype(string) ai(integer 0) ///
	Neighbor(integer 1) population altvariance exog(varlist fv) ///
	pscorr(integer 0) psfixnote(integer 0) dp(varname) xvars(varlist) ///
	selfxvars(varlist) vgamma(name)]

if (`pscorr') unab _n1list : _n1-_n`neighbor'

// return model info
if ("`exog'"!="") return local exog = "`exog'"

if ("`varlist'"!="") {
	return local depvar = "`varlist'"
}
else {
	exit
}

// create header output table
di as text "{hline 28}{c TT}{hline 59}"
di as text "        Variable     Sample {c |}    Treated     Controls   Difference         S.E.   T-stat"
di as text "{hline 28}{c +}{hline 59}"

// r(table) setup
tempname RTAB C
local coleq ""
local colnm ""
local _level = c(level)
local _crit  = invnormal(1 - (100 - `_level') / 200)

// create body and return results
qui foreach v of varlist `varlist' {
	// no matched outcome for obs off support
	replace _`v' = . if _support==0
	cap replace _self_`v' = . if _support==0

	tempname m1t m0t u0u u1u m0u m1u att atu seatt seatu seate

	sum `v' if _treated==1, mean
	scalar `u1u' = r(mean)
	sum `v' if _treated==0, mean
	scalar `u0u' = r(mean)

	sum `v' if _treated==1 & _support==1, mean
	scalar `m1t' = r(mean)
	local N1 = r(N)
	sum _`v' if _treated==1 & _support==1, mean
	scalar `m0t' = r(mean)
	scalar `att' = `m1t' - `m0t'

	if ("`ate'"!="") {
		sum _`v' if _treated==0 & _support==1, mean
		scalar `m1u' = r(mean)
		sum `v' if _treated==0 & _support==1, mean
		scalar `m0u' = r(mean)
		local N0 = r(N)
		scalar `atu' = `m1u' - `m0u'
		scalar `ate' = `att'*`N1'/(`N0'+`N1') + `atu'*`N0'/(`N0'+`N1')
	}

	if (`ai' != 0) {
		tempvar VhatE VhatEt VhatEu shat w
		g `w' = max(_weight, 0)
		// AI (2006, eq14 p. 250), or Aetal (2004, p303)
		g `shat' = cond("`altvariance'" == "", (`ai' / (`ai' + 1)) * (`v' - _self_`v')^2, _self_`v')
		if ("`population'" == "") { // AI (2006) Theorem 6 p.250: conditional variance
			g `VhatEt' = `shat' * (_treated - (1 - _treated) * `w')^2
			if ("`ate'" != "") {
				g `VhatEu' = `shat' * ((1 - _treated) - _treated * `w')^2
				g `VhatE'  = `shat' * (1 + `w')^2
			}
		}
		else {  // AI (2006) Theorem 7 p.251: marginal variance
			g `VhatEt' = max(0, _treated * (`v' - _`v' - `att')^2) ///
				+ (1 - _treated) * (`w'^2 - `w' / `neighbor') * `shat'
			if ("`ate'" != "") {
				g `VhatEu' = max(0, (1 - _treated) * (_`v' - `v' - `atu')^2) ///
					+ _treated * (`w'^2 - `w' / `neighbor') * `shat'
				g `VhatE'  = max(0, (_treated * (`v' - _`v') + (1 - _treated) * (_`v' - `v') - `ate')^2) ///
					+ (`w'^2 + 2 * `w' - `w' / `neighbor') * `shat'
				}
		}
		sum `VhatEt' if _support==1, mean
		scalar `seatt' = sqrt(r(sum)) / `N1'
		if ("`ate'"!="") {
			sum `VhatEu' if _support==1, mean
			scalar `seatu' = sqrt(r(sum)) / `N0'
			sum `VhatE'  if _support==1, mean
			scalar `seate' = sqrt(r(sum)) / (`N0' + `N1')
		}
	}
	else { // calculate approx s.e.'s
		tempname wtot var1 var0 number
		tempvar w2
		
		sum `v' if _treated==1 & _support==1  
		scalar `var1'   = r(Var)
		sum `v' if _treated==0 & _weight<.
		scalar `var0'   = r(Var)
		gen `w2' = _weight^2 if _treated==0
		sum `w2'
		scalar `wtot' = r(sum)
		if ("`spline'"!="") | ("`llr'"!="" & "`kerneltype'"=="tricube")   {
			scalar `seatt' = .
		}
		else scalar `seatt' = sqrt(`var1'/`N1' + `var0'*`wtot'/`N1'^2)
		scalar `seatu' = .
		scalar `seate' = .
	}

	if (`pscorr') {
		tempname seate_ai seatt_ai seatu_ai vv _psr
		scalar `seate_ai' = `seate'
		scalar `seatt_ai' = `seatt'
		scalar `seatu_ai' = `seatu'
		local _att_v = `att'
		local _atu_v = `atu'
		local _ate_v = `ate'
		mata: st_matrix("`_psr'", pscorr_ai2016("`v'", "_self_`v'", "_`v'", ///
			"`xvars'", "`selfxvars'", "`dp'", "`vgamma'", "`_n1list'", ///
			`ai', `neighbor', `N0', `N1', `_att_v', `_atu_v', `_ate_v'))
		// ATE
		scalar `vv' = `seate_ai'^2 - `_psr'[1,1]
		if (`vv' >= -1e-10 & `vv' < 0) scalar `vv' = 0
		if (`vv' >= 0) {
			scalar `seate' = sqrt(`vv')
		}
		else {
			noi di as text "Warning: corrected ATE variance < -1e-10 for `v'; SE set to missing."
			scalar `seate' = .
		}
		// ATT
		scalar `vv' = `seatt_ai'^2 - `_psr'[1,2] + `_psr'[1,3]
		if (`vv' >= -1e-10 & `vv' < 0) scalar `vv' = 0
		if (`vv' >= 0) {
			scalar `seatt' = sqrt(`vv')
		}
		else {
			noi di as text "Warning: corrected ATT variance < -1e-10 for `v'; SE set to missing."
			scalar `seatt' = .
		}
		// ATU
		scalar `vv' = `seatu_ai'^2 - `_psr'[1,4] + `_psr'[1,5]
		if (`vv' >= -1e-10 & `vv' < 0) scalar `vv' = 0
		if (`vv' >= 0) {
			scalar `seatu' = sqrt(`vv')
		}
		else {
			noi di as text "Warning: corrected ATU variance < -1e-10 for `v'; SE set to missing."
			scalar `seatu' = .
		}
	}

	// build r(table) columns for this outcome
	matrix `C' = J(9, 1, .)
	matrix `C'[1,1] = `att'
	matrix `C'[2,1] = `seatt'
	if (`seatt' < . & `seatt' > 0) {
		local _z = `att' / `seatt'
		matrix `C'[3,1] = `_z'
		matrix `C'[4,1] = 2 * normal(-abs(`_z'))
		matrix `C'[5,1] = `att' - `_crit' * `seatt'
		matrix `C'[6,1] = `att' + `_crit' * `seatt'
	}
	matrix `C'[8,1] = `_crit'
	matrix `C'[9,1] = 0
	matrix `RTAB' = nullmat(`RTAB'), `C'
	local coleq "`coleq' `v'"
	local colnm "`colnm' ATT"
	if ("`ate'" != "") {
		matrix `C' = J(9, 1, .)
		matrix `C'[1,1] = `atu'
		matrix `C'[2,1] = `seatu'
		if (`seatu' < . & `seatu' > 0) {
			local _z = `atu' / `seatu'
			matrix `C'[3,1] = `_z'
			matrix `C'[4,1] = 2 * normal(-abs(`_z'))
			matrix `C'[5,1] = `atu' - `_crit' * `seatu'
			matrix `C'[6,1] = `atu' + `_crit' * `seatu'
		}
		matrix `C'[8,1] = `_crit'
		matrix `C'[9,1] = 0
		matrix `RTAB' = nullmat(`RTAB'), `C'
		local coleq "`coleq' `v'"
		local colnm "`colnm' ATU"
		matrix `C' = J(9, 1, .)
		matrix `C'[1,1] = `ate'
		matrix `C'[2,1] = `seate'
		if (`seate' < . & `seate' > 0) {
			local _z = `ate' / `seate'
			matrix `C'[3,1] = `_z'
			matrix `C'[4,1] = 2 * normal(-abs(`_z'))
			matrix `C'[5,1] = `ate' - `_crit' * `seate'
			matrix `C'[6,1] = `ate' + `_crit' * `seate'
		}
		matrix `C'[8,1] = `_crit'
		matrix `C'[9,1] = 0
		matrix `RTAB' = nullmat(`RTAB'), `C'
		local coleq "`coleq' `v'"
		local colnm "`colnm' ATE"
	}

	tempname ols seols
	qui regress `v' _treated
	scalar `ols' = _b[_treated]
	scalar `seols' = _se[_treated]

	noi di as text %16s abbrev("`v'",16) "  Unmatched {c |}" as result %11.0g `u1u' "  " %11.0g `u0u' "  " %11.0g `ols' "  " %11.0g `seols' "  " %7.2f `ols'/`seols'
	noi di as text              _col(17) "        ATT {c |}" as result %11.0g `m1t' "  " %11.0g `m0t' "  " %11.0g `att'	"  " %11.0g `seatt' "  " %7.2f `att'/`seatt'

	if ("`ate'"!="") {
		noi di as text _col(17) "        ATU {c |}" as result %11.0g `m0u' "  " %11.0g `m1u' "  " %11.0g `atu'	"  " %11.0g `seatu' "  " %7.2f `atu'/`seatu'
		noi di as text _col(17) "        ATE {c |}" _col(56) as result %11.0g `ate'	"  " %11.0g `seate' "  " %7.2f `ate'/`seate' 
	}
	noi di as text "{hline 28}{c +}{hline 59}"

	// return estimates and s.e.'s
	return scalar att = `att'
	return scalar att_`v' = `att'
	return scalar seatt = `seatt'
	return scalar seatt_`v' = `seatt'
	if ("`ate'"!="") {
		return scalar atu = `atu'
		return scalar atu_`v' = `atu'
		return scalar ate = `ate'
		return scalar ate_`v' = `ate'
		return scalar seatu = `seatu'
		return scalar seatu_`v' = `seatu'
		return scalar seate = `seate'
		return scalar seate_`v' = `seate'
		if (`pscorr') {
			return scalar seate_ai_fixed_`v' = `seate_ai'
			return scalar seatt_ai_fixed_`v' = `seatt_ai'
			return scalar seatu_ai_fixed_`v' = `seatu_ai'
			return scalar qA_`v'      = `_psr'[1,1]
			return scalar qTminus_`v' = `_psr'[1,2]
			return scalar qTplus_`v'  = `_psr'[1,3]
			return scalar qUminus_`v' = `_psr'[1,4]
			return scalar qUplus_`v'  = `_psr'[1,5]
		}
	}

}

matrix rownames `RTAB' = b se z pvalue ll ul df crit eform
matrix coleq    `RTAB' = `coleq'
matrix colnames `RTAB' = `colnm'
return matrix table = `RTAB'

if (`ai'==0) {
	if (`seatt' != .) di as text "Note: S.E. does not take into account that the propensity score is estimated."
}
else if (`pscorr') {
	di as text "Note: Population S.E. adjusted for estimated propensity scores."
}
else {
	if ("`population'"=="") di as text "Note: Sample S.E."
	else di as text "Note: Population S.E."
	if (`psfixnote') di as text "Note: S.E. treats the propensity score as fixed."
}

tab _treated _support

end


// KERNEL MATCHING
program define _Match_kernel
	syntax anything [, OUTcome(varlist) Kerneltype(string) BWidth(real 0.06) CALiper(string) METric(string) W(string) ATE]
	tempname weight dif base
	tempvar out

	if ("`metric'"=="mahalanobis") {
		qui mata : _Dif_mbase("`anything'", "`base'", "`w'") 
	}

	count if _treated==0 & _support==1
	local N0 = r(N)
	if ("`ate'"!="") local start 1
	else local start = `N0' + 1

	count if _treated<=1 & _support==1
	forvalues obs = `start'/`r(N)' {
		_Dif_`metric' `anything' if _support==1 & _treated==(`obs'<=`N0'), obs(`obs') dif(`dif') base(`base') w(`w')
		_Kernel_ `kerneltype' `weight' `dif' `bwidth'
		replace _weight = _weight + `weight' if `weight'!=.
	
		if ("`outcome'"!="") {
			foreach v of varlist `outcome' {
				sum `v' [aw=`weight'] if _support==1 & _treated==(`obs'<=`N0'), mean
				if (r(mean)!=.) replace _`v' = r(mean) in `obs'
			}
		}
		cap assert `weight'==.
		if (_rc==0) replace _support = 0 in `obs'

		drop `weight' `dif'
	}
end


// LLR MATCHING
program define _Match_llr
	syntax anything [, OUTcome(varlist) Kerneltype(string) BWidth(real 0.06) CALiper(string) METric(string) W(string) ATE]
	tempname weight dif V base
	tempvar out
	
	if ("`metric'"=="mahalanobis") {
		mata : _Dif_mbase("`anything'", "`base'", "`w'")
	}

	count if _treated==0 & _support==1
	local N0 = r(N)
	if ("`ate'"!="") local start 1
	else local start = `N0' + 1

	count if _treated<=1 & _support==1
	forvalues obs = `start'/`r(N)' {
		_Dif_`metric' `anything' if _support==1 & _treated==(`obs'<=`N0'), obs(`obs') dif(`dif') base(`base') w(`w')
		_Kernel_ `kerneltype' `weight' `dif' `bwidth'
		sum `dif' [aw=`weight']
		scalar `V' = r(Var) * (r(N) - 1) / r(N)
		replace _weight = _weight + `weight' * (`V' + r(sum)^2 - r(sum) * `dif') / `V' if `weight'!=.
		if ("`outcome'"!="") {
			foreach v of varlist `outcome' {
				cap reg `v' `dif' [aw=`weight'] if _support==1 & _treated==(`obs'<=`N0')
				if (!_rc) replace _`v' = _b[_cons] in `obs'
				else replace _support = 0 in `obs'
			}
		}
		drop `weight' `dif'
	}
end


// COMMON SUPPORT FUNCTIONS
program define _Support_
	syntax varname [, level(real 100) untreated ate]
	if (`level'==100) {
		sum `varlist' if _treated==0, mean
		replace _support = 0 if (`varlist'<r(min) | `varlist'>r(max)) & _treated==1
		if ("`ate'"!="") {
			sum `varlist' if _treated==1, mean
			replace _support = 0 if (`varlist'<r(min) | `varlist'>r(max)) & _treated==0
		}
	}
	else {
		_Support_trim `varlist', level(`level')
		if ("`ate'"!="") {
			_Support_trim `varlist', level(`level') treated(0)
		}
	}
end


program define _Support_trim
	syntax varname [, level(real 100) treated(integer 1)]
	tempvar x0 y0
	kdensity `varlist' if _treated==(1-`treated'), nograph at(`varlist') gen(`x0' `y0')
	replace _support = 0 if `y0'==0 & _treated==`treated'
	if (`level'>0) {
		_pctile `y0' if _treated==`treated', p(`level')
		replace _support = 0 if `y0'<r(r1) & _treated==`treated'
	}
end


// DIFFERENCING FUNCTIONS
program define _Dif_pscore
	syntax varname [if], obs(int) dif(string) [base(string) W(string)]
	qui g double `dif' = (`varlist' - `varlist'[`obs']) `if'
end


program define _Dif_mahalanobis
	syntax anything [if], obs(int) dif(string) base(string) W(string)
	tempname b
	mata : st_matrix("`b'", st_data(`obs', tokens("`anything'")))  // row vector x[i] into b
	matrix colnames `b' = `anything'
	matrix `b' = `b'*`w'
	matrix score double `dif' = `b' `if'                          // x' * W * x[i]
	replace `dif' = `base' - 2 * `dif' + `base'[`obs'] `if'       // (x - x[i]) * W * (x - x[i])
end


// VARIOUS KERNELS
program define _Kernel_
	args kernel weight dif bwidth
	
	if ("`kernel'"=="epan") {
		qui g double `weight' = 1 - (`dif'/`bwidth')^2 if abs(`dif')<=`bwidth'
	}
	else if ("`kernel'"=="normal") {
		qui g double `weight' = normalden(`dif'/`bwidth')
	}
	else if ("`kernel'"=="biweight") {
		qui g double `weight' = (1 - (`dif'/`bwidth')^2)^2 if abs(`dif')<=`bwidth'
	}
	else if ("`kernel'"=="uniform") {
		qui g double `weight' = 1 if abs(`dif')<=`bwidth'
	}
	else if ("`kernel'"=="tricube") {
		qui g double `weight' = (1-abs(`dif'/`bwidth')^3)^3 if abs(`dif')<=`bwidth'
	}
	// normalize sum of weights to 1
	sum `weight', mean
	replace `weight' = `weight'/r(sum)
end


mata:

// calculates x'Wx used by mahalanobis metric, needs to be done only once
real scalar _Dif_mbase(string xvars, string base, string wmatrix)
{
	real matrix W, X
	real scalar i, j
	j = st_addvar("double", base)
	st_view(X=., ., tokens(xvars))
	W = st_matrix(wmatrix)
	for (i = 1; i<=st_nobs(); ++i) _st_store(i, j, X[i,.] * W * X[i,.]')
	return(j)
}


void match_mahalanobis(real scalar t0, real scalar t1, real scalar c0, real scalar c1, 
	real scalar neighbors, real scalar caliper, real scalar noreplace, real scalar ties, 
	string wmatrix, string xvars, string n1, string outvar, string moutvar, string ate)
{
	real scalar obs, nout, nmatch
	real matrix W, dist, smallest
	real vector id
	
	real scalar altvar
	altvar = (st_local("altvariance") != "")

	st_view(X=., ., tokens(xvars))
	st_view(N1=., ., tokens(n1))
	st_view(WEIGHT=., ., "_weight")
	st_view(NN=., ., "_nn")
	st_view(SUPPORT=., ., "_support")
	st_view(XWX=., ., _Dif_mbase(xvars, st_tempname(), wmatrix))
	nout = rows(tokens(outvar))
	if (nout > 0) {
		st_view(OUTVAR=., ., tokens(outvar)) 
		st_view(MOUTVAR=., ., tokens(moutvar))
	}

	id = st_data(.,"_id")   // we want a copy
	W = st_matrix(wmatrix)
	
	for (i = t0; i<=t1; ++i) {
		if (t0 == c0 && ate=="" && WEIGHT[i]==0) continue
		dist = XWX[c0..c1] - 2*X[c0..c1,.]*(W*X[i,.]') :+ XWX[i]
		dist = dist, id[c0..c1]
		if (t0 == c0) dist = select(dist, dist[.,2]:!=i)
		smallest = kth_smallest(dist, neighbors, caliper)
		nmatch = rows(smallest)
		if (t0 != c0) NN[i] = nmatch
		for (k = 1; k<=nmatch; ++k) {
			obs = smallest[k]
			if (t0 != c0) {
				N1[i, k] = obs
				WEIGHT[obs] = WEIGHT[obs] + 1 / neighbors
			}
			if (nout > 0) {
				MOUTVAR[i,.] = MOUTVAR[i,.] + OUTVAR[obs,.] :/ nmatch
			}
		}
		// estimate conditional variance following Abadie et al. (2004, p.303)
		if (altvar && t0 == c0) {
			m = (nmatch :* MOUTVAR[i,.] + OUTVAR[i,.]) :/ (nmatch + 1)
			MOUTVAR[i,.] = (OUTVAR[i,.] - m)^2
			for (k = 1; k <= nmatch; ++k) {
				MOUTVAR[i,.] = MOUTVAR[i,.] + (OUTVAR[smallest[k],.] - m)^2
			}
			MOUTVAR[i,.] = MOUTVAR[i,.] :/ nmatch
		}

		if (nmatch < 1 && t0 != c0) SUPPORT[i] = 0
	}

} // end of match_mahal


real matrix kth_smallest(real matrix a, real scalar k, real scalar caliper)
{
    real scalar i, j, l, m, x, tmp

    l = 1
	m = rows(a)
    while (l < m) {
      x = a[k, 1]
      i = l
      j = m
      do {
		while (a[i, 1] < x) i++
		while (x < a[j, 1]) j--
		if (i <= j) {
			tmp = a[i, .]
			a[i, .] = a[j, .]
			a[j, .] = tmp
			i++
			j--
		}
      } while (i <= j)
      if (j < k) l = i
      if (k < i) m = j
    }
	return(select(a[1..k,2], a[1..k,1] :<= caliper)) 
} // end of kth_smallest


void match_pscore(real scalar i0, real scalar i1, real scalar j0, real scalar j1, 
	real scalar neighbors, real scalar caliper, real scalar noreplace, real scalar ties, 
	string wmatrix, string xvars, string n1, string outvar, string moutvar, string ate)
{
	real scalar dif0, dif1, obs, i, jmatch, j, k, nmatch, nout, forward, idx_idlist, idx_ismatch, altvar

	st_view(PSCORE=., ., "_pscore")
	st_view(WEIGHT=., ., "_weight")
	st_view(TREATED=., ., "_treated")
	st_view(SUPPORT=., ., "_support")
	st_view(N1=., ., tokens(n1))
	st_view(NN=., ., "_nn")
	st_view(ID=., ., "_id")

	altvar = (st_local("altvariance") != "")

	idx_idlist = st_addvar("long", st_tempname())
	st_view(IDLIST=., ., idx_idlist)

	idx_ismatch = st_addvar("byte", st_tempname())
	st_view(ISMATCH=., ., idx_ismatch)
	for(i=1; i<=rows(ISMATCH); i++) ISMATCH[i] = 0

	nout = rows(tokens(outvar))
	if (nout>0) {
		st_view(OUTVAR=., ., tokens(outvar)) 
		st_view(MOUTVAR=., ., tokens(moutvar))
	}

	nmatch = 0
	forward = 1
	i = i0
	jmatch = j0
	while (i<=i1 && (jmatch>=j0 && jmatch<=j1)) {
		if (i==j0) ++jmatch
		if (i==jmatch) --jmatch
		j = jmatch
		dif1 = abs(PSCORE[i] - PSCORE[j])
		while (j>=j0 && j<j1) {
			j = next_unmatched(i, j, forward, noreplace, idx_ismatch)
			if (j<j0 || j>j1) break
			dif0 = dif1
			dif1 = abs(PSCORE[i] - PSCORE[j])
			if (dif1>dif0) j = j1
			if (dif1<dif0) jmatch = j
		}
		// update match and match-ID variables
		if (abs(PSCORE[i] - PSCORE[jmatch]) < caliper) {
			ISMATCH[jmatch] = 1
			++nmatch
			IDLIST[nmatch] = ID[jmatch]
			if (ties>0) {       // match remaining ties
				nmatch = match_ties(i, jmatch, j0, j1, nmatch, forward, noreplace, idx_idlist, idx_ismatch)
			}
			if (neighbors>1) {  // match remaining neighbors (1-to-many)
				nmatch = match2(i, jmatch, j0, j1, neighbors, caliper, nmatch, idx_idlist)
			}
			for(k = 1; k<=nmatch; k++) {
				obs = IDLIST[k]
				if (i0 != j0 & k<=neighbors) { // note that with ties we only keep k<=neighbors id's
					N1[i, k] = obs
				}
				if (i0 != j0) WEIGHT[obs] = WEIGHT[obs] + 1/nmatch
				if (nout>0) MOUTVAR[i,.] = MOUTVAR[i,.] + OUTVAR[obs,.]:/nmatch
			}
			if (i0 != j0) NN[i] = nmatch
			// estimate conditional variance following Abadie et al. (2004, p.303)
			if (altvar && i0 == j0) {
				m = (nmatch :* MOUTVAR[i,.] + OUTVAR[i,.]) :/ (nmatch + 1)
				MOUTVAR[i,.] = (OUTVAR[i,.] - m)^2
				for (k = 1; k <= nmatch; ++k) {
					MOUTVAR[i,.] = MOUTVAR[i,.] + (OUTVAR[IDLIST[k],.] - m)^2
				}
				MOUTVAR[i,.] = MOUTVAR[i,.] :/ nmatch
			}


		} else if (i0 != j0) SUPPORT[i] = 0
	
		if (jmatch==j1 && i<i1) forward = 0
		if (noreplace==1 && nmatch>0) {
			jmatch = next_unmatched(i, jmatch, forward, noreplace, idx_ismatch)
		}
		nmatch = 0
		++i
	}
}


real scalar next_unmatched(real scalar obs, real scalar j0, real scalar forward, real scalar noreplace, real scalar idx_ismatch)
{
	real scalar j
	j = j0
	do {
		if (forward==1) ++j
		else --j
		
		if (obs!=j & noreplace==0) return(j)
	} while (obs==j | _st_data(j, idx_ismatch)==1)
	return(j)
}


real scalar match2(real scalar obs, real scalar jm, real scalar j0, real scalar j1, real scalar neighbors, real scalar caliper, real scalar nmatch, real scalar idx_idlist)
{
	real scalar dif0, dif1, k, pos0, pos1, jmatch

	st_view(PSCORE=., ., "_pscore")
	st_view(ID=., ., "_id")
	st_view(IDLIST=., ., idx_idlist)

	jmatch = jm
	k = 1
	pos0 = jmatch - 1
	if (pos0==obs) --pos0
	pos1 = jmatch + 1
	if (pos1==obs) ++pos1
	while (k<neighbors) {
		if (pos0>=j0) dif0 = abs(PSCORE[obs] - PSCORE[pos0])
		if (pos1<=j1) dif1 = abs(PSCORE[obs] - PSCORE[pos1])
		if ((dif0<=dif1 && pos0>=j0) || (pos0>=j0 && pos1>j1)) {
			jmatch = pos0
			--pos0
			if (pos0==obs) --pos0
		} else if (pos1<=j1) {
			jmatch = pos1
			++pos1
			if (pos1==obs) ++pos1
		} else k = neighbors
		// update match and match-ID variables
		if (abs(PSCORE[obs] - PSCORE[jmatch])<caliper && k<neighbors) {
			++nmatch
			IDLIST[nmatch] = ID[jmatch]
			++k
		} else k = neighbors
	}
	return(nmatch)
}


real scalar match_ties(real scalar obs, real scalar jmatch, real scalar j0, real scalar j1, real scalar nmatch, real scalar forward, real scalar noreplace, real scalar idx_idlist, real scalar idx_ismatch)
{
	real scalar i, dif0, dif

	st_view(PSCORE=., ., "_pscore")
	st_view(ID=., ., "_id")
	st_view(IDLIST=., ., idx_idlist)
	st_view(ISMATCH=., ., idx_ismatch)

	dif0 = abs(PSCORE[obs] - PSCORE[jmatch])

	i = next_unmatched(obs, jmatch, forward, noreplace, idx_ismatch)
	if (i<j0 || i>j1) return(nmatch)
	
	dif = abs(PSCORE[obs] - PSCORE[i])
	while (dif<=dif0) {
		++nmatch
		IDLIST[nmatch] = ID[i]
		if (noreplace==1) ISMATCH[i] = 1

		i = next_unmatched(obs, i, forward, noreplace, idx_ismatch)
		if (i<j0 || i>j1) return(nmatch)
		
		dif = abs(PSCORE[obs] - PSCORE[i])
	}
	return(nmatch)
}


// Abadie-Imbens (2016) first-stage correction for propensity-score matching.
//
// For ATT/ATU, the published variance formula contains two score-estimation
// terms.  The first is the c_t (or c_u) term in Theorem 2.  The second is the
// derivative of the target parameter with respect to the propensity-score
// parameter, d tau_t(theta)/d theta.
//
// Abadie and Imbens (2016, p. 799) estimate this derivative by matching on the
// full covariate vector X:
//
//   E[p(X)]^{-1} E[ X f(X'theta) { mu(1,X) - mu(0,X) - tau_t } ].
//
// This implementation uses the equivalent conditional-on-propensity-score
// population decomposition
//
//   E[p(X)]^{-1} E[ f { E[X|p] (mubar_1(p)-mubar_0(p)-tau_t)
//      + cov(X,mu(1,X)|p) - cov(X,mu(0,X)|p) } ].
//
// Thus dT_vec and dU_vec are alternative plug-in estimators of the same
// population derivative, not literal implementations of the full-X matching
// estimator displayed on p. 799.  The same-arm local means used below are
// leave-one-out through selfY/selfX.
real rowvector pscorr_ai2016(
	string scalar yvar,      string scalar selfy_var, string scalar matchy_var,
	string scalar xvars_str, string scalar selfxvars_str,
	string scalar dp_var,    string scalar vgamma_name,
	string scalar n1vars_str,
	real scalar ai_val,      real scalar M,
	real scalar N0,          real scalar N1,
	real scalar att_val,     real scalar atu_val,   real scalar ate_val)
{
	real matrix Y, selfY, matchY, X, selfX, dP, PS, TR, SP, N1mat, VG
	string matrix cstripe
	real scalar n_obs, K, K_x, cons_pos, factor, i, k, m, mid, n_matches
	real scalar fi, pi, qi, tri, mu1_i, mu0_i, aT_i, aU_i
	real matrix Cown, Cother, C1, C0, Cfull1, Cfull0, Xfull
	real matrix cA, cT, dT_vec, cU, dU_vec
	real scalar qA, qTminus, qTplus, qUminus, qUplus

	st_view(Y=., .,      yvar)
	st_view(selfY=., .,  selfy_var)
	st_view(matchY=., ., matchy_var)
	st_view(X=., .,      tokens(xvars_str))
	st_view(selfX=., .,  tokens(selfxvars_str))
	st_view(dP=., .,     dp_var)
	st_view(PS=., .,     "_pscore")
	st_view(TR=., .,     "_treated")
	st_view(SP=., .,     "_support")
	st_view(N1mat=., .,  tokens(n1vars_str))

	VG    = st_matrix(vgamma_name)
	K     = cols(VG)
	K_x   = cols(tokens(xvars_str))
	n_obs = rows(Y)

	// locate _cons column in Vgamma
	cstripe  = st_matrixcolstripe(vgamma_name)
	cons_pos = 0
	for (k = 1; k <= K; k++) {
		if (cstripe[k, 2] == "_cons") {
			cons_pos = k
			break
		}
	}

	if (K != K_x + (cons_pos > 0)) {
		_error(3498, "pscorr_ai2016: X and Vgamma dimensions do not match")
	}

	// Cown_i = (ai/(ai+1)) * (X_i - selfX_i) * (Y_i - selfY_i)
	factor = ai_val / (ai_val + 1)
	Cown   = J(n_obs, K_x, 0)
	for (i = 1; i <= n_obs; i++) {
		if (SP[i,1] != 1) continue
		Cown[i,.] = factor :* (X[i,.] - selfX[i,.]) :* (Y[i,1] - selfY[i,1])
	}

	cA     = J(1, K, 0)
	cT     = J(1, K, 0)
	dT_vec = J(1, K, 0)
	cU     = J(1, K, 0)
	dU_vec = J(1, K, 0)

	for (i = 1; i <= n_obs; i++) {
		if (SP[i,1] != 1) continue

		fi  = dP[i,1]
		pi  = PS[i,1]
		qi  = 1 - pi
		if (pi <= 0 | pi >= 1) {
			_error(3498, "pscorr_ai2016: propensity score outside (0,1)")
		}
		tri = TR[i,1]

		// average Cown over cross-arm matches (_n1..._nM)
		n_matches = 0
		Cother = J(1, K_x, 0)
		for (m = 1; m <= M; m++) {
			mid = N1mat[i, m]
			if (mid >= .) continue
			n_matches++
			Cother = Cother + Cown[mid,.]
		}
		if (n_matches == 0) _error(3498, "pscorr_ai2016: no match IDs for obs " + strofreal(i))
		Cother = Cother / n_matches

		if (tri == 1) {
			C1 = Cown[i,.]
			C0 = Cother
		}
		else {
			C0 = Cown[i,.]
			C1 = Cother
		}

		// insert _cons slot: 0 in Cfull, 1 in Xfull
		if (cons_pos == 0) {
			Cfull1 = C1
			Cfull0 = C0
			Xfull  = X[i,.]
		}
		else if (cons_pos == 1) {
			Cfull1 = (0, C1)
			Cfull0 = (0, C0)
			Xfull  = (1, X[i,.])
		}
		else if (cons_pos == K) {
			Cfull1 = (C1, 0)
			Cfull0 = (C0, 0)
			Xfull  = (X[i,.], 1)
		}
		else {
			Cfull1 = (C1[1,1..cons_pos-1], 0, C1[1,cons_pos..K_x])
			Cfull0 = (C0[1,1..cons_pos-1], 0, C0[1,cons_pos..K_x])
			Xfull  = (X[i,1..cons_pos-1],  1, X[i,cons_pos..K_x])
		}

		// ATE
		cA = cA + fi :* (Cfull1 :/ pi + Cfull0 :/ qi)

		// Local propensity-score means for ATT/ATU. selfY is leave-one-out.
		if (tri == 1) {
			mu1_i = selfY[i,1]
			mu0_i = matchY[i,1]
		}
		else {
			mu0_i = selfY[i,1]
			mu1_i = matchY[i,1]
		}

		aT_i = mu1_i - mu0_i - att_val
		aU_i = mu1_i - mu0_i - atu_val

		// cT and cU implement the Abadie-Imbens (2016) c_t/c_u terms.
		// dT_vec and dU_vec implement the derivative terms using the
		// conditional-on-propensity decomposition documented above:
		//   X*a + C1 - C0.
		// This avoids a separate full-X matching pass for the derivative
		// component but targets the same population derivative.
		cT     = cT     + fi :* (Xfull :* aT_i + Cfull1 + (pi/qi) :* Cfull0)
		dT_vec = dT_vec + fi :* (Xfull :* aT_i + Cfull1 - Cfull0)
		cU     = cU     + fi :* (Xfull :* aU_i - Cfull0 - (qi/pi) :* Cfull1)
		dU_vec = dU_vec + fi :* (Xfull :* aU_i + Cfull1 - Cfull0)
	}

	cA     = cA     :/ (N0 + N1)
	cT     = cT     :/ N1
	dT_vec = dT_vec :/ N1
	cU     = cU     :/ N0
	dU_vec = dU_vec :/ N0

	qA      = (cA     * VG * cA')[1,1]
	qTminus = (cT     * VG * cT')[1,1]
	qTplus  = (dT_vec * VG * dT_vec')[1,1]
	qUminus = (cU     * VG * cU')[1,1]
	qUplus  = (dU_vec * VG * dU_vec')[1,1]

	return((qA, qTminus, qTplus, qUminus, qUplus))
}


end
