# EpiStrainDynamics (development version)

## Minor improvements

* Plot functions now call `ggplot2` functions explicitly rather than importing
  the whole package (#35).
* Standardised quotation-mark style (#37).
* Reordered function arguments so those without defaults precede those with
  defaults (#38).
* Enabled parallel testing with `Config/testthat/parallel: true` (#39).
* Added tests for MCMC validation.
* Standardised error handling on `cli::cli_abort()` (#49).

## Documentation

* Replaced `\dontrun{}` in examples with `@examplesIf interactive()` for slow
  examples and guards such as `@examplesIf rlang::is_installed("xts")` where a
  package is used only in an example (#36).
* Corrected documentation of both `sarscov2` and `influenza` datasets (#56).
* Clarify contribution and maintenance guidelines, installation from source
  instructions, and add citations to readme (#51, #57, #58, #59, #60).
* Expanded documentation for `smoothing_structure()` and `diagnose_model()`, 
  add cross-linked references where applicable (#53, #54, #55).
* Increased R version dependency in line with base pipe usage (#43).

# EpiStrainDynamics 0.0.1 (2026-03-27)

Initial version submitted to rOpenSci for peer review
(ropensci/software-review#763).

## New features

* Random-walk and p-spline models for inferring temporal trends of single and
  multiple (including subtyped) pathogens, fitted via Stan.
* Model construction and fitting via `construct_model()` and `fit_model()`,
  with `smoothing_structure()` and `dispersion_structure()` helpers.
* Post-processing generics `incidence()`, `growth_rate()`, `Rt()`, and
  `proportion()`, each with a `plot()` method.
* `diagnose_model()` for MCMC convergence diagnostics.
* Bundled `influenza` and `sarscov2` example datasets.
* "Using EpiStrainDynamics" and "Algorithmic Scaling" vignettes.
