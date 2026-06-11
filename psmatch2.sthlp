{smcl}
{hline}
help for {hi:psmatch2}
{hline}

{title:Mahalanobis and Propensity score Matching}

{p 8 21 2}{cmdab:psmatch2}
{it:depvar}
[{it:indepvars}]
[{cmd:if} {it:exp}]
[{cmd:in} {it:range}]
[{cmd:,}
    {cmdab:out:come}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:p:score}{cmd:(}{it:varname}{cmd:)}
    {cmdab:n:eighbor}{cmd:(}{it:integer}{cmd:)}
    {cmd:radius}
    {cmdab:cal:iper}{cmd:(}{it:real}{cmd:)}
    {cmdab:mahal:anobis}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:ai}{cmd:(}{it:integer}{cmd:)}
    {cmd:samplevar}
    {cmdab:altv:ariance}
    {cmd:kernel}
    {cmd:llr}
    {cmdab:k:erneltype}{cmd:(}{it:type}{cmd:)}
    {cmdab:bw:idth}{cmd:(}{it:real}{cmd:)}
    {cmd:spline}
    {cmdab:n:knots}{cmd:(}{it:integer}{cmd:)}
    {cmdab:com:mon}
    {cmd:trim}{cmd:(}{it:real}{cmd:)}
    {cmdab:norepl:acement}
    {cmdab:desc:ending}
    {cmd:odds}
    {cmd:index}
    {cmd:logit}
    {cmd:ties}
    {cmdab:qui:etly}
    {cmd:w}{cmd:(}{it:matrix}{cmd:)}
    {cmd:ate}]
    
{p 8 21 2}where {it:indepvars} and {cmdab:mahal:anobis}{cmd:(}{it:varlist}{cmd:)} may contain factor variables; see {cmd:fvvarlist}.

{title:Description}

{pstd}
{cmd:psmatch2} implements Mahalanobis matching and several propensity-score
matching estimators. The command compares treated and untreated observations
after adjusting for pre-treatment observed covariates. Treatment status is
defined by {it:depvar}=1 for treated observations and {it:depvar}=0 for
untreated observations.

{pstd}
The propensity score is the conditional probability of treatment. It may be
supplied by the user with {cmd:pscore()} or estimated internally from
{it:indepvars} by probit or logit. Matching methods include nearest-neighbor
matching, k-nearest-neighbor matching, caliper matching, radius matching,
kernel matching, local linear regression matching, spline smoothing, and
Mahalanobis matching.

{pstd}
{cmd:psmatch2} is being continuously improved. To install the most recent SSC
version, type

{phang2}{cmd:. ssc install psmatch2, replace}{p_end}

{title:Standard errors}

{pstd}
By default, {cmd:psmatch2} reports approximate standard errors that treat the
matching weights as fixed. These standard errors assume independent
observations, homoskedastic outcome variances within the treated and untreated
groups, and outcome variances that do not vary with the propensity score. They
also do not account for estimation of the propensity score.

{pstd}
For the ATT, the default standard error is

{phang2}
SE(ATT) = sqrt( Var(Y|DM=1)/N1 + Var(Y|DM=0)*sum(w_i^2; i in DM=0)/N1^2 )
{p_end}

{pstd}
where N1 is the number of matched treated observations, DM=1 denotes the
matched treated sample, DM=0 denotes the matched controls, and w_i is the
number of times control observation i is used as a match. With option
{cmd:ate}, analogous formulas are used for the ATU and ATE:

{phang2}
SE(ATU) = sqrt( Var(Y|DM=0)/N0 + Var(Y|DM=1)*sum(w_i^2; i in DM=1)/N0^2 )
{p_end}

{phang2}
SE(ATE) = sqrt( [ Var(Y|DM=1)*sum((1+w_i)^2; i in DM=1) + Var(Y|DM=0)*sum((1+w_i)^2; i in DM=0) ] / (N0+N1)^2 )
{p_end}

{pstd}
where N0 is the number of matched controls and w_i is the number of times each
observation is used as a match by the opposite treatment arm.

{pstd}
For nearest-neighbor matching, analytical Abadie-Imbens (2006) standard errors
are available with {cmd:ai(}{it:M}{cmd:)}, where {it:M}>0. The value {it:M} is
the number of neighbors used to estimate the conditional outcome variance
σ²(X,W), as in their formula (14).

{pstd}
By default, {cmd:ai()} reports the marginal, or population, variance of the
matching estimator (Theorem 7 of Abadie and Imbens 2006). Specify
{cmd:samplevar} to report the conditional/sample variance instead (Theorem 6).
Under {cmd:samplevar}, the variance convention conditions on the realized
matching sample, and the estimated propensity score is treated as fixed.

{pstd}
When the propensity score is estimated internally by probit or logit and
{cmd:ai()} is specified, {cmd:psmatch2} applies the Abadie-Imbens (2016)
correction for first-step score estimation whenever the correction is
available. This correction adjusts the analytical AI standard errors for the
fact that the propensity score is estimated rather than known.

{pstd}
The AI(2016) correction is not applied when the propensity score is supplied
with {cmd:pscore()} or when {cmd:samplevar}, {cmd:caliper}, {cmd:ties},
{cmd:noreplacement}, or {cmd:altvariance} is specified.

{pstd}
Bootstrapping nearest-neighbor matching estimators is generally not
recommended. Use {cmd:ai(}{it:M}{cmd:)} for Abadie-Imbens analytical matching
standard errors.

{title:Returned results}

{pstd}
{cmd:psmatch2} stores the ATT in {cmd:r(att)} and its standard error in
{cmd:r(seatt)}. With option {cmd:ate}, it also stores the ATU and ATE in
{cmd:r(atu)} and {cmd:r(ate)}, with standard errors in {cmd:r(seatu)} and
{cmd:r(seate)}. With more than one outcome variable, outcome-specific results
are stored as {cmd:r(att_}{it:varname}{cmd:)},
{cmd:r(seatt_}{it:varname}{cmd:)}, and similarly for ATU and ATE.

{pstd}
{cmd:psmatch2} also stores {cmd:r(table)}, a Stata-style results table with
rows {cmd:b}, {cmd:se}, {cmd:z}, {cmd:pvalue}, {cmd:ll}, {cmd:ul}, {cmd:df},
{cmd:crit}, and {cmd:eform}. The columns identify the reported outcome-effect
combinations.

{pstd}
Option {cmd:ate} controls which treatment-effect parameters are reported. It
does not determine the standard error for the ATT. Thus, when the AI(2016)
correction is available, the ATT standard error is the same whether or not
{cmd:ate} is specified.

{title:Replication and ties}

{pstd}
To make results replicable, set the random-number seed before calling
{cmd:psmatch2}. The sort order of the data may affect nearest-neighbor matching
when there are ties in the propensity score, for example when the propensity
score is estimated from categorical covariates.

{pstd}
You can open the dialog by {dialog psmatch2:clicking here} or by typing

{phang2}{cmd:. db psmatch2}{p_end}

{pstd}
The following list presents the syntax for each matching method.

{title:Sample weights}

{pstd}
{cmd:psmatch2} estimates the propensity score in the unweighted analysis
sample and uses the estimated score to construct matches. Sample weights
can be used after matching to compute matched outcomes or matched treatment effects
for a target population.

{pstd}
For a population ATT, compare the weighted mean of the observed treated outcome
with the weighted mean of the matched counterfactual outcome among treated
observations on support:

{phang2}{cmd:. sum outcome if treated==1 & _support==1 [aw=pweight]}{p_end}
{phang2}{cmd:. sum _outcome if treated==1 & _support==1 [aw=pweight]}{p_end}

{pstd}
The difference between these two weighted means is the population-weighted ATT
computed from the matched outcomes left behind by {cmd:psmatch2}.

{pstd}
For a population ATU, use the sampling weights for the untreated population
and compare the weighted mean of {cmd:_outcome} with the weighted mean of
{cmd:outcome} among untreated observations on support.

{pstd}
For a population ATE, average the matched individual treatment effects over all
observations on support using the sampling weights. One direct calculation is:

{phang2}{cmd:. gen double _te = cond(_treated==1, outcome - _outcome, _outcome - outcome) if _support==1}{p_end}
{phang2}{cmd:. sum _te if _support==1 [aw=pweight]}{p_end}

{pstd}
For weighted estimands, assess covariate balance using the sampling weights.
If weighted balance is poor, change the propensity-score specification or the
matching rule, for example by adding sampling-design variables or interactions
to the propensity-score model.

{title:Matching within strata}

{pstd} The following code illustrates how to match within exact cells and then calculate the average effect for the whole population.

	{cmd:g att = .}
	{cmd:egen g = group(groupvars)}
	{cmd:levels g, local(gr)}
	{cmd:qui foreach j of local gr {c -(}}
		{cmd:psmatch2 treatvar varlist if g==`j', out(outvar)}
		{cmd:replace att = r(att) if  g==`j'}
	{cmd:{c )-}}
	{cmd:sum att}

{title:Detailed Syntax}

{phang}
{bf:One-to-one matching:}

{p 8 21 2}{cmdab:psmatch2}
{it:depvar}
[{it:indepvars}]
[{cmd:if} {it:exp}]
[{cmd:in} {it:range}]
{cmd:,}
    [{cmdab:out:come}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:p:score}{cmd:(}{it:varname}{cmd:)}
    {cmdab:ai}{cmd:(}{it:integer k}>1{cmd:)}
    {cmdab:mahal:anobis}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:cal:iper}{cmd:(}{it:real}{cmd:)}
    {cmdab:norep:lacement}
    {cmdab:desc:ending}
    {cmdab:com:mon}
    {cmd:trim}{cmd:(}{it:real}{cmd:)}
    {cmd:odds}
    {cmd:index}
    {cmd:logit}
    {cmd:ties}
    {cmdab:warn:ings}
    {cmdab:qui:etly}
    {cmd:ate}]


{phang}
{bf:{it:k}-Nearest neighbors matching:}

{p 8 21 2}{cmdab:psmatch2}
{it:depvar}
[{it:indepvars}]
[{cmd:if} {it:exp}]
[{cmd:in} {it:range}]
{cmd:,}
    [{cmdab:out:come}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:p:score}{cmd:(}{it:varname}{cmd:)}
    {cmdab:n:eighbor}{cmd:(}{it:integer k}>1{cmd:)}
    {cmdab:cal:iper}{cmd:(}{it:real}{cmd:)}
    {cmdab:com:mon}
    {cmd:trim}{cmd:(}{it:real}{cmd:)}
    {cmd:odds}
    {cmd:index}
    {cmd:logit}
    {cmd:ties}
    {cmdab:warn:ings}
    {cmdab:qui:etly}
    {cmd:ate}]


{pstd}
{bf:Radius matching:}

{p 8 21 2}{cmdab:psmatch2}
{it:depvar}
[{it:indepvars}]
[{cmd:if} {it:exp}]
[{cmd:in} {it:range}]
{cmd:,}
    {cmdab:radius}
    {cmdab:cal:iper}{cmd:(}{it:real}{cmd:)}
    [{cmdab:out:come}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:p:score}{cmd:(}{it:varname}{cmd:)}
    {cmdab:com:mon}
    {cmd:trim}{cmd:(}{it:real}{cmd:)}
    {cmd:odds}
    {cmd:index}
    {cmd:logit}
    {cmdab:qui:etly}
    {cmd:ate}]


{pstd}
{bf:Kernel matching:}

{p 8 21 2}{cmdab:psmatch2}
{it:depvar}
[{it:indepvars}]
[{cmd:if} {it:exp}]
[{cmd:in} {it:range}]
{cmd:,}
    {cmdab:kernel}
    [{cmdab:out:come}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:k:erneltype}{cmd:(}{it:kernel_type}{cmd:)}
    {cmdab:p:score}{cmd:(}{it:varname}{cmd:)}
    {cmdab:bw:idth}{cmd:(}{it:real}{cmd:)}
    {cmdab:mahal:anobis}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:com:mon}
    {cmd:trim}{cmd:(}{it:real}{cmd:)}
    {cmd:odds}
    {cmd:index}
    {cmd:logit}
    {cmdab:qui:etly}
    {cmd:ate}]


{pstd}
{bf:Local linear regression matching:}

{p 8 21 2}{cmdab:psmatch2}
{it:depvar}
[{it:indepvars}]
[{cmd:if} {it:exp}]
[{cmd:in} {it:range}]
{cmd:,}
    {cmdab:llr}
    {cmdab:out:come}{cmd:(}{it:varlist}{cmd:)}
    [{cmdab:k:erneltype}{cmd:(}{it:kernel_type}{cmd:)}
    {cmdab:p:score}{cmd:(}{it:varname}{cmd:)}
    {cmdab:bw:idth}{cmd:(}{it:real}{cmd:)}
    {cmdab:mahal:anobis}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:com:mon}
    {cmd:trim}{cmd:(}{it:real}{cmd:)}
    {cmd:odds}
    {cmd:index}
    {cmd:logit}
    {cmdab:qui:etly}
    {cmd:ate}]

{phang}
{bf:Spline matching:}

{p 8 21 2}{cmdab:psmatch2}
{it:depvar}
[{it:indepvars}]
[{cmd:if} {it:exp}]
[{cmd:in} {it:range}]
{cmd:,}
    {cmd:spline}
    {cmdab:out:come}{cmd:(}{it:varlist}{cmd:)}
    [{cmdab:nk:nots(}{it:integer}{cmd:)}
    {cmdab:p:score}{cmd:(}{it:varname}{cmd:)}
    {cmdab:n:eighbor}{cmd:(}{it:integer}{cmd:)}
    {cmdab:cal:iper}{cmd:(}{it:real}{cmd:)}
    {cmdab:com:mon}
    {cmd:trim}{cmd:(}{it:real}{cmd:)}
    {cmd:odds}
    {cmd:index}
    {cmd:logit}
    {cmd:ties}
    {cmdab:warn:ings}
    {cmdab:qui:etly}
    {cmd:ate}]

{pstd}
{bf:Mahalanobis matching:}

{p 8 21 2}{cmdab:psmatch2}
{it:depvar}
[{cmd:if} {it:exp}]
[{cmd:in} {it:range}]
{cmd:,}
    {cmdab:mahal:anobis}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:out:come}{cmd:(}{it:varlist}{cmd:)}
    {cmdab:ai}{cmd:(}{it:integer}{cmd:)}
    [{cmd:samplevar}
    {cmdab:altv:ariance}
    {cmdab:k:ernel}{cmd:(}{it:kernel_type}{cmd:)}
    {cmd:llr}
    {cmdab:bw:idth}{cmd:(}{it:real}{cmd:)}
    {cmdab:cal:iper}{cmd:(}{it:real}{cmd:)}
    {cmd:w}{cmd:(}{it:matrix}{cmd:)}
    {cmd:ate}]


{pstd}
{cmd:psmatch2} creates a number of variables for the convenience of the user:

{pmore}
{inp:_treated} is a variable that equals 0 for control observations and 1 for treatment observations.

{pmore}
{inp:_support} is an indicator variable with equals 1 if the observation is on the common support
and 0 if the observatio is off the support.

{pmore}
{inp:_pscore} is the estimated propensity score or a copy of the one provided by {cmdab:p:score()}.

{pmore}
{inp:_{it:outcome_variable}} for every treatment observation stores the value of the matched outcome.

{pmore}
{inp:_weight}. For nearest neighbor matching, it holds the frequency with which the
observation is used as a match; with option {cmd:ties} and k-nearest neighbors matching it holds
the normalized weight; for kernel matching, and llr matching with a weight other than
stata's tricube, it stores the overall weight given to the matched observation. When estimating att only
_weight = 1 for the treated.

{pmore}
{inp:_id} In the case of one-to-one and nearest-neighbors matching, a new identifier created for all observations.

{pmore}
{inp:_n{it:k}} In the case of one-to-one and nearest-neighbors matching, for every treatment observation,
it stores the observation number of the k-th matched control observation. Do not forget to sort by
{it:_id} if you want to use the observation number (id) of for example the 1st nearest neighbor as in

	{cmd:sort _id}
	{cmd:g x_of_match = x[_n1]}

{pmore}
{inp:_nn} In the case of nearest-neighbors matching, for every treatment observation,
it stores the number of matched control observations.

{title:Options}

{phang}
{cmdab:out:come}{cmd:(}{it:varlist}{cmd:)} the outcome variable(s). 
When evaluating multiple outcomes psmatch2 reduces to the min common number of observations 
with non-missing values on ALL outcomes, because otherwise the matching weigths will not sum 
to the right number. If you have multiple outcomes with widely differing missing values 
you may wish to run psmatch2 separately for each of the outcomes.

{phang}
{cmdab:ate} with this option the average treatment effect (ate) and average treatment
effect on the untreated (atu) are reported in addition to the average treatment effect on the treated (att).
The estimates are returned in {it:r(ate)}, {it:r(atu)} and {it:r(att)} respectively, see above.

{phang}
{cmdab:ai}{cmd:(}{it:integer}{cmd:)}
with nearest-neighbor matching, calculate the Abadie-Imbens analytical
standard errors
proposed by Abadie and Imbens (2006) by specifying the number of neighbors {it:M} used
to estimate the conditional outcome variance σ²(X,W) (their formula (14)). With option {cmdab:altv:ariance} one can
use the estimator of Abadie et al. (2004) instead.

{pmore}
By default, {cmd:ai()} reports the marginal, or population, variance of the
matching estimator (Theorem 7 of Abadie and Imbens 2006), which adds a
component for treatment-effect heterogeneity across covariate values. Specify
{cmd:samplevar} to report the conditional/sample variance instead (Theorem 6),
which conditions on the realized matching sample. Under {cmd:samplevar}, the
estimated propensity score is treated as fixed, so the AI(2016) first-stage
correction is not applied.

{pmore}
For Mahalanobis matching ({cmd:mahal()}), the AI(2006) standard errors are returned directly.

{pmore}
For propensity-score nearest-neighbor matching with an internally estimated
probit or logit score, {cmd:ai()} additionally applies the Abadie-Imbens
(2016) correction for first-stage score estimation when the correction is
available. The correction is not applied when the score is supplied with
{cmd:pscore()} or when
{cmd:caliper}, {cmd:ties}, {cmd:noreplacement}, {cmd:altvariance},
{cmd:common}, {cmd:index}, {cmd:odds}, {cmd:samplevar}, {cmd:kernel},
{cmd:llr}, {cmd:radius}, {cmd:spline}, or {cmd:mahalanobis} is specified.
The ATE correction is weakly negative in variance. For ATT and ATU, the correction can increase or decrease the SE.
The option {cmd:ate} determines whether ATU and ATE are displayed and returned.
It does not determine whether the ATT first-stage correction is applied.

{pmore}
Implementation note. For ATT and ATU, Abadie and Imbens (2016, p. 799) estimate
the derivative of the target parameter with respect to the propensity-score
parameter by matching on the full covariate vector. {cmd:psmatch2} instead uses
an equivalent population decomposition conditional on the propensity score. In
particular, the derivative term is written as the sum of a local
propensity-score mean component and the difference between the within-score
covariance terms for the treated and untreated outcome regressions. This avoids
a separate full-covariate matching step for this component. The same-arm local
means are computed leave-one-out. Thus the ATT and ATU first-stage corrections
are plug-in implementations of the Abadie-Imbens population correction, but the
derivative component is not the literal full-covariate matching estimator
displayed in their paper.

{phang}
{cmdab:samplevar}
when using {cmdab:ai}{cmd:(}{it:integer}{cmd:)}, report the conditional/sample
variance of the matching estimator (Theorem 6 of Abadie and Imbens 2006)
rather than the marginal, or population, variance (Theorem 7, the default).
The population variance adds a component for treatment-effect heterogeneity,
V^τ(X), on top of the conditional variance. Under {cmd:samplevar}, the
estimated propensity score is treated as fixed and the AI(2016) correction is
not applied.

{phang}
{cmdab:altv:ariance}
when using {cmdab:ai}{cmd:(}{it:integer}{cmd:)}, calculate the conditional variance using the expression in Abadie et al. (2004, p.303).

{title:Options: Estimation of the propensity score}

{phang}
{cmdab:p:score}{cmd:(}{it:varname}{cmd:)} specifies the variable to be used as propensity score.

{pstd}
Alternatively, {it:indepvars} need to be specified to allow the program to estimate the propensity score on them.
In this case:

{phang}
{cmd:logit} use logit instead of the default probit to estimate the propensity score.

{phang}
{cmdab:qui:etly} do not print output of propensity score estimation.

{phang}
{cmd:odds} match on the logarithm of the odds ratio of the propensity score (stored in _pscore).

{phang}
{cmd:index} use the latent variable index instead of the probability.

{phang}
{cmdab:warn:ings} test for control observations with duplicate propensity score values.


{title:Options: Imposition of common support}

{phang}
{cmdab:com:mon} imposes a common support by dropping treatment observations
whose pscore is higher than the maximum or less than the minimum pscore of the controls.

{phang}
{cmd:trim(}{it:integer}{cmdab:)} imposes common support by dropping # percent of the treatment
observations at which the pscore density of the control observations is the lowest.

{title:Options: Choice of matching estimator}

{phang}
{cmdab:n:eighbor}{cmd:(}{it:integer}{cmd:)} number of neighbors used to calculate the matched outcome.
Defaults to 1. Default matching method is single nearest-neighbour (without caliper).

{phang}
{cmdab:norep:lacement} perform 1-to-1 matching without replacement. Nearest neigbor propensity score matching only.

{phang}
{cmdab:desc:ending} perform 1-to-1 matching without replacement in descending order. Nearest neighbor propensity score matching only.

{phang}
{cmd:ties} not only match nearest neighbor but also other controls with identical (tied) pscores.

{phang}
{cmd:radius} perform radius matching within the specified radius given by {cmd:caliper}.

{phang}
{cmdab:cal:iper}{cmd:(}{it:real}{cmd:)} value for maximum distance of controls.
Use to perform nearest neighbor(s) within caliper, radius matching and Mahalanobis 1-to-1 matching.

{phang}
{cmdab:k:ernel} perform kernel matching.

{phang}
{cmdab:k:erneltype}{cmd:(}{it:kernel_type}{cmd:)} specifies the type of kernel:

{p 8 8 2}{inp:normal} the gaussian kernel.

{p 8 8 2}{inp:biweight} the biweight kernel.

{p 8 8 2}{inp:epan} the epanechnikov kernel (Default).

{p 8 8 2}{inp:uniform} the uniform kernel.

{p 8 8 2}{inp:tricube} the tricube kernel.

{phang}
{cmd:llr} llr use local linear regression matching instead of kernel matching.

{phang}
{cmdab:bw:idth}{cmd:(}{it:real}{cmd:)} the bandwidth for kernel and local linear regression matching.
Default bandwidth is 0.06, except when doing local linear regression with the Epanechnikov kernel when the default bandwidth is the rule-of-thumb bandwidth of {cmd:lpoly}.

{phang}
{cmdab:mahal:anobis}{cmd:(}{it:varlist}{cmd:)} perform Mahalanobis-metric matching on {it: varlist}.

{phang}
{cmd:w}{cmd:(}{it:matrix}{cmd:)} specify alternative weighting matrix. Mahalanobis-metric matching
becomes matching on a quadratic metric with the specified weighting matrix.

{phang}
{cmd:spline} performs 'spline-smoothing matching' by first fitting a natural cubic
    spline on pscore (or on the result from estimate) to outcome.
    The matched values are stored in the new variable, _s_outcomevar.  (It requires the
    {cmd:spline} programme, which for stata7 needs to be downloaded by typing: net install snp7_1.)

{phang}
{cmdab:nk:nots(}{it:integer}{cmd:)} specifies the number of interior knots for spline smoothing. Default is
    the fourth root of the number of comparison units.

{title:Saved results}

{pstd}
{cmd:psmatch2} saves the following in {cmd:r()}:

{synoptset 28 tabbed}{...}
{synopt:{cmd:r(att)}}average treatment effect on the treated{p_end}
{synopt:{cmd:r(seatt)}}standard error of ATT{p_end}
{synopt:{cmd:r(ate)}}average treatment effect (with {cmd:ate}){p_end}
{synopt:{cmd:r(seate)}}standard error of ATE (with {cmd:ate}){p_end}
{synopt:{cmd:r(atu)}}average treatment effect on the untreated (with {cmd:ate}){p_end}
{synopt:{cmd:r(seatu)}}standard error of ATU (with {cmd:ate}){p_end}
{synopt:{cmd:r(table)}}Stata-style results table with rows {cmd:b}, {cmd:se}, {cmd:z}, {cmd:pvalue}, {cmd:ll}, {cmd:ul}, {cmd:df}, {cmd:crit}, and {cmd:eform}{p_end}

{pstd}
With multiple outcome variables, effects and standard errors are also returned
as outcome-specific scalars such as {cmd:r(att_}{it:varname}{cmd:)},
{cmd:r(seatt_}{it:varname}{cmd:)}, {cmd:r(ate_}{it:varname}{cmd:)}, and
{cmd:r(seate_}{it:varname}{cmd:)}.

{pstd}
When the Abadie-Imbens (2016) first-stage correction fires, the following
additional scalars are returned for each outcome variable {it:y}:

{synoptset 28 tabbed}{...}
{synopt:{cmd:r(seatt_ai_fixed_}{it:y}{cmd:)}}AI(2006) SE for ATT before the AI(2016) correction{p_end}
{synopt:{cmd:r(qTminus_}{it:y}{cmd:)}, {cmd:r(qTplus_}{it:y}{cmd:)}}ATT correction terms: first-stage covariance term and derivative term{p_end}
{synopt:{cmd:r(seate_ai_fixed_}{it:y}{cmd:)}}AI(2006) SE for ATE before the correction (with {cmd:ate}){p_end}
{synopt:{cmd:r(seatu_ai_fixed_}{it:y}{cmd:)}}AI(2006) SE for ATU before the correction (with {cmd:ate}){p_end}
{synopt:{cmd:r(qA_}{it:y}{cmd:)}}correction term for ATE (with {cmd:ate}): {cmd:r(seate)}^2 = {cmd:r(seate_ai_fixed_}{it:y}{cmd:)}^2 - {cmd:r(qA_}{it:y}{cmd:)}{p_end}
{synopt:{cmd:r(qUminus_}{it:y}{cmd:)}, {cmd:r(qUplus_}{it:y}{cmd:)}}ATU correction terms: first-stage covariance term and derivative term (with {cmd:ate}){p_end}

{title:Examples}

    {inp: . psmatch2 training age gender, kernel k(biweight) out(wage)}
    {inp: . psmatch2 training age gender, n(5) logit}
    {inp: . psmatch2 training age gender, out(wage)}
    {inp: . bs "psmatch2 training age gender, out(wage)" "r(att)"}

{title:Also see}

{pstd}
The commands {help pstest}, {help psgraph}.


{title:Thanks for citing {cmd:psmatch2} as follows}

{pstd}
E. Leuven and B. Sianesi. (2003). "PSMATCH2: Stata module to perform full Mahalanobis and propensity score matching, common support graphing, and covariate imbalance testing".
http://ideas.repec.org/c/boc/bocode/s432001.html. This version INSERT_VERSION_HERE.

where you can check your version as follows:

    {inp: . which psmatch2}


{title:Disclaimer}

{pstd}
THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

{pstd}
IN NO EVENT WILL THE COPYRIGHT HOLDERS OR THEIR EMPLOYERS, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR REDISTRIBUTE THIS SOFTWARE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

{title:Background Reading}

{p 0 2}Abadie, A., Drukker, D., Herr, J. L., & Imbens, G. W. (2004). "Implementing matching estimators for average treatment effects in Stata", {it:Stata journal 4}, 290-311.

{p 0 2}Abadie, A. and Imbens, G.W. (2006), "Large sample properties of matching estimators for average treatment effects", {it:Econometrica 74}(1), 235-267.

{p 0 2}Abadie, A. and Imbens, G.W. (2016), "Matching on the Estimated Propensity Score", {it:Econometrica 84}(2), 781-807.

{p 0 2}Cochran, W. and Rubin, D.B. (1973), "Controlling Bias in Observational Studies", {it:Sankyha 35}, 417-446.

{p 0 2}Dehejia, R.H and Wahba, S. (1999), "Causal Effects in Non-Experimental Studies: Re-Evaluating the Evaluation of Training Programmes", {it:Journal of the American Statistical Association 94}, 1053-1062.

{p 0 2}Heckman, J.J., Ichimura, H. and Todd, P.E. (1997), "Matching As An Econometric Evaluation Estimator: Evidence from Evaluating a Job Training Programme", {it:Review of Economic Studies 64}, 605-654.

{p 0 2}Heckman, J.J., Ichimura, H. and Todd, P.E. (1998), "Matching as an Econometric Evaluation Estimator", {it:Review of Economic Studies 65}, 261-294.

{p 0 2}Heckman, J.J., Ichimura, H., Smith, J.A. and Todd, P. (1998), "Characterising Selection Bias Using Experimental Data", {it:Econometrica 66}, 5.

{p 0 2}Heckman, J.J., LaLonde, R.J., Smith, J.A. (1998), "The Economics and Econometrics of Active Labour Market Programmes", in Ashenfelter, O. and Card, D. (eds.), {it:The Handbook of Labour Economics Vol. 3A}.

{p 0 2}Imbens, G. (2000), "The Role of Propensity Score in Estimating Dose-Response Functions", {it:Biometrika 87(3)}, 706-710.

{p 0 2}Lechner, M. (2001), Identification and Estimation of Causal Effects of Multiple Treatments under the Conditional Independence Assumption, in: Lechner, M., Pfeiffer, F. (eds), {it:Econometric Evaluation of Labour Market Policies}, Heidelberg: Physica/Springer, p. 43-58.

{p 0 2}Rosenbaum, P.R. and Rubin, D.B. (1983), "The Central Role of the Propensity Score in Observational Studies for Causal Effects", {it:Biometrika 70}, 1, 41-55.

{p 0 2}Rosenbaum, P.R. and Rubin, D.B. (1985), "Constructing a Control Group Using Multivariate Matched Sampling Methods that Incorporate the Propensity Score", {it:The American Statistician 39(1)}, 33-38.

{p 0 2}Rubin, D.B. (1974), "Estimating Causal Effects of Treatments in Randomised and Non-Randomised Studies", {it:Journal of Educational Psychology 66}, 688-701.

{p 0 2}Rubin, D.B. (1980), "Bias Reduction Using Mahalanobis-Metric Matching", {it:Biometrics 36}, 293-298.

{title:Author}

{pstd}
Edwin Leuven, University of Oslo. If you observe any problems {browse "mailto:e.leuven@gmail.com"}.

{pstd}
Barbara Sianesi, Institute for Fiscal Studies, London, UK.
