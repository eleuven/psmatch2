*! psmatch2_collect.ado  v1.0  2026-06-10
*  Posts r(table) from the most recent psmatch2 run into collect.
*  Call psmatch2 first; this program does not call psmatch2 itself.

program define psmatch2_collect
	version 17
	syntax [, CLEAR]

	if ("`clear'" != "") collect clear

	tempname T
	matrix `T' = r(table)

	local eqs : coleq `T'
	local cn  : colnames `T'
	local K = colsof(`T')

	forvalues j = 1/`K' {
		local y   : word `j' of `eqs'
		local eff : word `j' of `cn'

		collect get ///
			b      = (`T'[1,`j']) ///
			se     = (`T'[2,`j']) ///
			z      = (`T'[3,`j']) ///
			pvalue = (`T'[4,`j']) ///
			ll     = (`T'[5,`j']) ///
			ul     = (`T'[6,`j']), ///
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
