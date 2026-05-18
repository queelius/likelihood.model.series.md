## R CMD check results

0 errors | 0 warnings | 1 note

* One NOTE: "unable to verify current time" (system clock check, benign).

## Update (v0.9.3 -> v0.10.0)

This release reshapes the exp_series functional API (dexp_series,
pexp_series, qexp_series, rexp_series, hazard_exp_series,
surv.exp_series, mean.exp_series) into thin wrappers over
dist.structure::exp_series(), the canonical constructor for series-of-
exponential distributions. The aggregate-rate math now ships from one
source of truth across the rlang ecosystem.

The mean.exp_series method remains polymorphic so that legacy
structure(rates, class = "exp_series") code paths continue to work
alongside the new dist.structure list shape. A small .onLoad hook
re-registers our S3 method after dist.structure's namespace is loaded,
ensuring the polymorphic dispatch wins consistently.

No public API removals. The functional API remains exported for
backwards compatibility; existing user code is unaffected.

## Coordinated submission

This is part of a coordinated 3-package submission. All packages are
maintained by me. Updated versions being submitted simultaneously:

- flexhaz 0.5.2
- serieshaz 0.2.0
- maskedcauses 0.10.0 (this package)

## Test environments

* local Ubuntu 24.04, R 4.3.3
* win-builder (R-devel and R-release)
