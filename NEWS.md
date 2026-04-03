# McMiso 0.2.0

## Breaking Changes

* The default prior hyperparameters for `PBclassifier()` and `miso()` have
  changed from `a0 = 0.25, b0 = 0.25` to `a0 = 0.5, b0 = 0.5` (Jeffreys
  prior). Results from previous versions can be reproduced by explicitly
  passing `a0 = 0.25, b0 = 0.25`.

* `mcmiso()` is now an exported function. Previously it was internal. Its
  interface has also changed: the `future` parallel plan must now be set by
  the caller via `future::plan()` before calling `mcmiso()`, rather than
  being set internally. This follows best practices for the `future` package.

* The return value of `miso()` now has class `"miso"`, enabling a dedicated
  `print` method. Code that tested `class(fit) == "list"` will need updating.

## New Features

* `misoN()`: multivariable isotonic regression for **continuous** outcomes
  using a Normal-Inverse-Chi-Squared conjugate model.

* `mcmisoN()`: parallel computing wrapper for `misoN()` (requires the
  `future` package and a parallel plan set by the caller).

* `mcPBclassifier()`: parallel computing wrapper for `PBclassifier()`.

* `boundary()`: extracts the decision boundary (minimal positive set) from
  a fitted `"pbc"` object.

* New `print` methods: `print.pbc()`, `print.miso()`, `print.misoN()`,
  `print.boundary()`.

* Comprehensive input validation has been added to all exported functions.

## Performance Improvements

* `getScenesV3()` has been rewritten from a brute-force `expand.grid(2^K)`
  enumeration to a recursive backtracking algorithm. This substantially
  reduces memory usage and computation time for large numbers of unique
  feature combinations.

* Vectorized computation replaces inner loops in `SweepCombTogBinom()`,
  `SweepMcCombBinom()`, `SweepCombTogNorm()`, and `SweepMcCombNorm()`.

## Dependency Changes

* `dplyr` has been removed as a dependency (no longer used internally).

* `future` has been moved from `Imports` to `Suggests`. It is only required
  for the parallel computing wrappers (`mcmiso`, `mcPBclassifier`,
  `mcmisoN`). These functions check for the package at run time and give an
  informative error if it is not installed.

# McMiso 0.1.2

* Initial CRAN release.
