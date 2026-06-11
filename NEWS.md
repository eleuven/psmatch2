# News

## v4.1.0 (2026-06-11)

### New features

- **Abadie-Imbens (2016) first-stage correction.** When the propensity score
  is estimated internally by probit or logit and nearest-neighbor matching is
  used without caliper, ties, noreplacement, altvariance, or common-support
  trimming, `psmatch2` automatically applies the Abadie-Imbens (2016)
  correction that adjusts the analytical AI standard errors for first-stage
  estimation of the propensity score. The correction fires for the ATT
  regardless of whether the `ate` option is specified. Use `ai(#)` to request
  AI standard errors; the correction is applied whenever eligible.

- **Factor variables in the propensity-score model.** The `ai()` standard
  errors and the AI(2016) first-stage correction now support factor-variable
  specifications such as `i.x` and `i.x#c.z` in the propensity-score model.

- **Approximate ATU and ATE standard errors.** When `ai()` is not used or the
  AI(2016) correction is not available, `psmatch2` now computes approximate
  standard errors for the ATU and ATE using formulas analogous to the ATT
  approximation. These were previously returned as missing.

- **`samplevar` option.** Requests the conditional (sample) variance of the
  matching estimator (Theorem 6 of Abadie and Imbens 2006) instead of the
  marginal (population) variance (Theorem 7, the default). Under `samplevar`,
  the estimated propensity score is treated as fixed and the AI(2016)
  correction is not applied.

- **`debug` option.** Returns diagnostic components of the AI(2016) correction:
  the pre-correction AI(2006) standard errors and the correction terms
  (`qA`, `qTminus`, `qTplus`, `qUminus`, `qUplus`). These satisfy the variance
  identities documented in the help file. Intended for testing and
  reproducibility checks.

- **`r(table)`.** `psmatch2` now returns a Stata-style results matrix with rows
  `b`, `se`, `z`, `pvalue`, `ll`, `ul`, `df`, `crit`, and `eform`. Columns
  identify outcome-effect combinations and are compatible with `collect` and
  `etable`.

- **`pstest` is now r-class.** `pstest` stores its results in `r()`, including
  balance statistics and counts.

### Changes

- The marginal (population) variance is now the default with `ai()`. The
  previous default was the conditional (sample) variance. Use `samplevar` to
  restore the old behavior.

- The `population` option has been removed. Population variance is now always
  the default and does not need to be requested explicitly.

- The `ai()` option now requires nearest-neighbor matching. Specifying `ai()`
  with kernel, LLR, radius, or spline matching produces an error.

### Help file

- New sections document the SE formulas for the default approximation and for
  AI(2006), the conditions for the AI(2016) correction, returned results, and
  weight-based post-matching calculations.
