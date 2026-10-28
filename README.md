<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://gcore.jsdelivr.net/gh/ChengHu2024/cade@main/branding/logo-for-dark.png">
    <img src="https://gcore.jsdelivr.net/gh/ChengHu2024/cade@main/branding/logo-for-light.png" width="400" alt="cade">
  </picture>
</p>
# cade

<!-- badges: start -->
[![R-CMD-check](https://github.com/ChengHu2024/cade/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChengHu2024/cade/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT) <!-- License metadata note: see NEWS.md -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23609650.svg)](https://doi.org/10.5281/zenodo.23609650)
[![Version](https://img.shields.io/badge/version-6.0.0-informational.svg)](https://doi.org/10.5281/zenodo.23609650)
<!-- badges: end -->

**Composition Attribution of Differential Expression**

Separate an observed bulk transcriptomic difference into what a change in
cell-type composition would produce, and what changed *within* the cells.

```r
observed log2FC  =  composition term  +  within-cell residual
```

The composition term is computed from a single-cell reference held constant
across samples:

```
composition(g) = log2( Σ_k f_case[k] · mu[k,g] / Σ_k f_ctrl[k] · mu[k,g] )
```

The residual is a **bookkeeping identity, not a causal claim**. It equals the
within-cell effect only if two conditions hold: the reference matches the
tissue, and deconvolution is unbiased. CADE exists to make those two
conditions checkable rather than assumed.

---

## Install

```r
# from GitHub
devtools::install_github("ChengHu2024/cade")

# or from the archived tarball
install.packages("cade_6.0.0.tar.gz", repos = NULL, type = "source")

library(cade)
self_test()
```

`self_test()` takes about a second, needs no data, and checks four things:

```
[PASS] composition term recovers the truth        measured 1.0000  (> 0.99)
[PASS] injected within-cell change recovered      measured 0.9984  (in (0.90, 1.10))
[PASS] no false positives above tolerance         measured 0.0000  (< 0.02)
[PASS] reporting does not alter the decomposition measured 1.0000  (== 1, bit-identical)
```

The fourth criterion exists because the report layer must never touch a number.
If it passes, the implementation does what this document says.

---

## Quick start

```r
library(cade)

## 1. Check the installation
self_test()

## 2. Prepare four things
##    obs     named numeric: log2 fold change per gene
##    frac    samples x cell types: estimated proportions
##    mu      cell types x genes: per-cell means, LINEAR scale
##    is_case logical vector over samples
diagnose_inputs(obs, frac, mu, is_case)

## 3. Decompose
d <- decompose_composition(obs, frac, mu, is_case)

## 4. Ask whether the answer can be interpreted
check_calibration(d) |> calibration_verdict()

## 5. Attribute and report
report(d, file = "cade_report.md")
```

---

## The four checks that matter

CADE refuses to run, or warns, in each of these situations. All four
correspond to failure modes that are **silent** in a naive implementation.

### 1. Is the reference linear?

`mu` must be a per-cell mean on a linear scale. A log-scaled matrix passed by
mistake produces a meaningless ratio and, without a check, no error at all.
`check_inputs()` compares `max(mu)` against the range a CP10k reference
occupies and refuses if it looks log-scaled.

### 2. Is the calibration informative?

Regressing `obs` on the composition term gives a slope and an R². Both are
reported at three expression strata, and the verdict is stated in words. Two
rules are enforced:

- always stratify by **expression quantile**, never by `|logFC|` — filtering on
  the outcome variable makes the criterion circular;
- the verdict on the **top-25% stratum** is the one that counts, because
  low-expression genes dilute the slope.

### 3. Is the residual interpretable?

`attribute_genes()` assigns one of five classes. The one that matters most is
`over-predicted`: a residual whose sign is **opposite** to `obs` means the
composition term overshot, and in that case `pct = resid/obs` is not a
meaningful share. A naive implementation reports a large `pct` here and
misattributes the gene to composition.

### 4. Do the cell-type columns line up?

`frac` and `mu` must be aligned by cell-type **name**, checked with
`setequal()` **before** any subsetting. Checking after `f <- f[, kinds]` is a
tautology and catches nothing. This exact bug occurred during the development
of this package; it is documented in `inst/CHANGELOG-2026-10-24.md`.

---

## What it looks like

![portrait](https://gcore.jsdelivr.net/gh/ChengHu2024/cade@main/branding/portrait.png)

A real, unedited console session — nothing in these images is hand-written:

| | |
|---|---|
| ![session-A](https://gcore.jsdelivr.net/gh/ChengHu2024/cade@main/branding/session-A.png) | ![session-B](https://gcore.jsdelivr.net/gh/ChengHu2024/cade@main/branding/session-B.png) |
| The normal path | The package refusing bad input |

---

## Two references cannot be compared by output alone

`compare_references()` reports how far two references' results diverge.
`recommend_reference()` then gives guidance — but deliberately **does not pick a
winner from the data**.

A reference that produces a smoother calibration is not thereby the correct
one; choosing on that basis is how the choice becomes circular. The function
evaluates only properties that can be checked *without* looking at the output:
whether the reference contains normal tissue, and whether any cell type is
starved. If both references pass those checks and still disagree, its verdict
is that **the disagreement is the finding** — report both, do not select.

---

## Known limits

These are properties of the method, not bugs. They are stated here because a
tool that hides its failure modes is worse than no tool.

- **A reference contaminated by 10–20% goes undetected.** In a known-answer
  simulation, mixing 20% of each cell type's expression into the others caused
  57.6% of unchanged genes to be misclassified, while the calibration verdict
  still read "informative".
- **The false-positive rate for low-expression genes depends on how the
  observation is computed.** The same simulated data gave 19% false positives
  with `log2(bulk + 1)` and a much lower rate with an exact `log2` of the mean.
- **A within-cell change in a cell-type-specific gene is diluted.** The same
  change is recovered with a magnitude that varies by up to 11-fold depending
  on where the gene is expressed.
- **Systematic deconvolution bias is absorbed into the residual.** A 20%
  systematic overestimate of one case-group fraction lowered recovery to 0.778
  while calibration still read "informative".

---

## Repository layout

```
R/            the package source (10 files)
man/          documentation (21 pages)
tests/        143 unit tests
inst/         cheat sheet (Chinese) and changelog
dev/          the scripts that produced the figures and the verification
branding/     figures used in this README
```

`dev/` is included on purpose: `verify_authenticity.R` is the five-part check
that the report layer does not alter any number, and `stress_test.R` produces
the known-answer results quoted above. You can run both yourself.

---

## Citation

```bibtex
@software{cade2026,
  author    = {Hu, Cheng},
  title     = {CADE: composition attribution of differential expression
               (R package, v6.0.0)},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23609650},
  url       = {https://doi.org/10.5281/zenodo.23609650}
}
```

The method is described in a preprint:

> Hu C. The identifiability boundary of composition attribution in bulk
> transcriptomics: a methodological study. Zenodo; 2026.
> doi:10.5281/zenodo.23655967

The first application of CADE - a four-cohort composition-attribution analysis
in abdominal aortic aneurysm - is deposited with its frozen data package:

> Hu C. Frozen analysis package for the identifiability boundary of composition
> attribution in bulk transcriptomics. Zenodo; 2026.
> doi:10.5281/zenodo.23656096

A companion manuscript on six failure modes is included in the same Zenodo
record as the methodology paper.

---

## Licence

MIT © 2026 Cheng Hu

If you find a case where this package is wrong, please tell me:
2024184006@hrbmu.edu.cn · ORCID [0009-0009-4627-0229](https://orcid.org/0009-0009-4627-0229)
