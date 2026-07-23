# EpiStrainDynamics (development version)

## Bug fixes

* Fixed the negative-binomial likelihood in the random-walk single-pathogen
  model when `dow_effect = TRUE` (#42).
* Corrected the `proportion()` / `plot()` example in the vignette and fixed
  the underlying figure path so plots render correctly (#45).
* Fixed the `Rt()` equation failing to render on the pkgdown site by setting
  MathJax as the math renderer (#44).

## Minor improvements

* Plot functions now call `ggplot2` functions explicitly rather than importing
  the whole package (#35).
* Standardised quotation-mark style (#37).
* Reordered function arguments so those without defaults precede those with
  defaults (#38).
* Enabled parallel testing with `Config/testthat/parallel: true` (#39).
* Added tests for MCMC validation.
* Added a `print()` method for `EpiStrainDynamics.model` objects (#47).
* Standardised error handling on `cli::cli_abort()` (#49).
* `plot()` is no longer redefined as a generic (#48).
* Increased the minimum R version to 4.1.0, matching actual base pipe usage,
  and added it to the R-CMD-check.yaml test matrix (#43).
* Fixed all lintr-identified issues (#50): styling (line length, indentation,
  spacing) via `styler`, dropped explicit `return()` on final statements,
  converted `sapply()` to `vapply()` for type-safe extraction, and removed
  genuine dead code flagged by unused-variable checks.

## Documentation

* Replaced `\dontrun{}` in examples with `@examplesIf interactive()` for slow
  examples and guards such as `@examplesIf rlang::is_installed("xts")` where a
  package is used only in an example (#36).
* Expanded documentation for `Rt()`, including its `gi_dist` argument (#44).
* Expanded documentation for `smoothing_structure()` and `diagnose_model()`
  (#53, #54).
* Expanded documentation for `growth_rate()`, `incidence()`, and
  `proportion()` to match the level of detail in `Rt()`, including symbol
  breakdowns for their adjustment formulas and, for `growth_rate()`, its
  relationship to `Rt()` (#66).
* Added roxygen `[]` cross-links for referenced functions throughout the
  documentation (#55).
* Committed `precompile.R` for the vignette and included revision
  instructions in CONTRIBUTING.md (#46).
* Clarified contribution and maintenance guidelines, added installation-
  from-source instructions, and added citations to the README (#51, #57,
  #58, #59, #60).
* Corrected documentation of both `sarscov2` and `influenza` datasets (#56).
* Added explanations of "pathogen noise" and the day-of-week effect (#61), a
  link to the relevant timeseries class (#62), and a COVID-19 example (#63),
  to the vignette.
* Clarified the G1.6 srr standard statement to note that, although
  performance/scaling is evaluated in the vignettes, this does not
  constitute an explicit comparison against alternative implementations
  (#52).

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
