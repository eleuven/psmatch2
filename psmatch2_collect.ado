*! psmatch2_collect.ado  v1.2  2026-06-11
*  Posts r(table) from the most recent psmatch2 run into collect.
*  Call psmatch2 first; this program does not call psmatch2 itself.

// Helper: loads one column of r(table) into r() scalars for collect get.
// Called before each collect get call so that r(b) r(se) etc. are available.
capture program drop _psmatch2_post_row
program define _psmatch2_post_row, rclass
	args mat j
	return scalar b      = `mat'[1,`j']
	return scalar se     = `mat'[2,`j']
	return scalar z      = `mat'[3,`j']
	return scalar pvalue = `mat'[4,`j']
	return scalar ll     = `mat'[5,`j']
	return scalar ul     = `mat'[6,`j']
end

program define psmatch2_collect
	version 17
	syntax [, CLEAR]

	capture confirm matrix r(table)
	if (_rc) {
		di as err "r(table) not found; run psmatch2 before psmatch2_collect"
		exit 301
	}

	// Copy r(table) before any rclass command can clear r()
	tempname T
	matrix `T' = r(table)

	// A repeated call otherwise leaves duplicate outcome#effect#result
	// items in different cmdsets, making the layout ambiguous.
	collect clear

	local eqs : coleq `T'
	local cn  : colnames `T'
	local K = colsof(`T')

	forvalues j = 1/`K' {
		local y   : word `j' of `eqs'
		local eff : word `j' of `cn'

		// Post this column into r() so collect get can read it
		_psmatch2_post_row `T' `j'
		collect get r(b) r(se) r(z) r(pvalue) r(ll) r(ul), ///
			tags(outcome[`y'] effect[`eff'])
	}

	collect label dim outcome "Outcome", modify
	collect label dim effect  "Effect",  modify

	collect label levels result ///
		b      "Coef."     ///
		se     "Std. err." ///
		z      "z"         ///
		pvalue "P>|z|"     ///
		ll     "Lower"     ///
		ul     "Upper",    modify

	collect layout (outcome#effect) (result[b se z pvalue ll ul])
end
