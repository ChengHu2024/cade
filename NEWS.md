# cade 6.0.0

Released 2026-10-25.

## Breaking changes

- **`compose_prediction()` and `decompose_composition()` gain a `convention`
  argument, and the default is now `"mean_of_log"`.** The previous behaviour
  (`"log_of_mean"`) is still available but must be requested explicitly.
  These are not interchangeable: on GSE232911 the two conventions give
  comp-obs correlations of 0.2846 and 0.4251, a difference of 0.14. Code that
  relied on the old default will produce different `pred` values. The choice
  is recorded on the result as `attr(x, "convention")`.

- **`cor_ci()` is new and has no defaults for `resample` or `unit`.** A caller
  who does not state what was resampled gets an error rather than an
  interval. The two units give intervals that differ by a factor of three to
  four in practice.

## New functions

- **`null_permutation()`** — exact permutation test of the comp-obs
  correlation under a null that scrambles the cell-type-to-profile
  correspondence. Reports `r_obs`, the full null distribution, `z`,
  `n_ge`, `p` and its Monte Carlo standard error. Warns when `n_perm` is too
  small to stabilise a `p` near 0.05. Cross-validated against an independent
  NumPy implementation on GSE232911: `r_obs` agrees to 4.8e-07 and `p`
  agrees within Monte Carlo error.

- **`cor_ci()`** — sampling uncertainty of the comp-obs correlation with an
  explicit resampling unit. `resample = "proportions"` re-estimates
  proportions only; `"proportions_and_obs"` also recomputes the observed log
  fold changes from per-sample expression and therefore requires `expr`. The
  returned object carries a `caveat` string that states in words what the
  interval does and does not cover.

- **`calibrate_by_layer()`** — fits the `obs ~ pred` calibration separately
  within each anatomical layer, with expression strata computed inside the
  layer. Reports a verdict per layer; a layer whose slope interval covers zero
  is marked `not_calibrated` rather than dropped.

- **`verify_deposit()`** — deposit-integrity checks A1 (intensity vs
  Detection P rank correlation per column, using `|rho|` so the check is
  indifferent to the sign convention), A2 (platform fingerprint against a
  reference deposit) and A3 (additive identity `S = log2(I) + r_i + c_j`).
  When a reference is supplied it also attempts to solve the row permutation
  and reports `permutation_recoverable`. The function states numerical
  verdicts only and makes no inference about cause.

- **`granularity_crossed()`** — runs the reference-granularity crossed design
  over true levels x reference levels x aggregation weights x seeds. Marks a
  reference level `not identifiable` when the spread of `r` across conditions
  exceeds a threshold. **Its `r` is computed without deconvolution**, so
  absolute values are not comparable with the companion paper; it tests the
  direction and magnitude of how the spread grows as the reference is
  coarsened.

## Bug fixes

- **`mean_of_log` no longer masks zeros with `pmax(..., 1e-12)`.** Masking
  turned an undefined ratio into a large finite number and silently defeated
  the zero-denominator guard that `log_of_mean` has. Both conventions now
  return `NA` for undefined genes.

- Column alignment in `compose_prediction()` asserts set equality *before*
  subsetting; the earlier form asserted after subsetting and could not detect
  a reference mismatch.

## Tests

- 120 test blocks at version 5.3.0; 170 at 6.0.0, all passing.
- `test-nullperm.R`, `test-deposit.R`, `test-granularity.R` are new.
- `test-mutation-killers.R` now exercises both `convention` values
  explicitly and asserts that the two do not coincide on a fixture whose
  within-group rows differ.

## Documentation

- `roxygen2` regenerated `NAMESPACE`; five exports and four `print` methods
  added.
- Every new function documents what it does **not** do, in a dedicated
  paragraph.
