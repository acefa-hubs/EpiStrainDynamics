# Changelog

## EpiStrainDynamics (development version)

### Bug fixes

- Fixed the negative-binomial likelihood in the random-walk
  single-pathogen model when `dow_effect = TRUE`
  ([\#42](https://github.com/acefa-hubs/EpiStrainDynamics/issues/42)).
- Corrected the
  [`proportion()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/proportion.md)
  /
  [`plot()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/plot.md)
  example in the vignette and fixed the underlying figure path so plots
  render correctly
  ([\#45](https://github.com/acefa-hubs/EpiStrainDynamics/issues/45)).
- Fixed the
  [`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md)
  equation failing to render on the pkgdown site by setting MathJax as
  the math renderer
  ([\#44](https://github.com/acefa-hubs/EpiStrainDynamics/issues/44)).

### Minor improvements

- Plot functions now call `ggplot2` functions explicitly rather than
  importing the whole package
  ([\#35](https://github.com/acefa-hubs/EpiStrainDynamics/issues/35)).
- Standardised quotation-mark style
  ([\#37](https://github.com/acefa-hubs/EpiStrainDynamics/issues/37)).
- Reordered function arguments so those without defaults precede those
  with defaults
  ([\#38](https://github.com/acefa-hubs/EpiStrainDynamics/issues/38)).
- Enabled parallel testing with `Config/testthat/parallel: true`
  ([\#39](https://github.com/acefa-hubs/EpiStrainDynamics/issues/39)).
- Added tests for MCMC validation.
- Added a [`print()`](https://rdrr.io/r/base/print.html) method for
  `EpiStrainDynamics.model` objects
  ([\#47](https://github.com/acefa-hubs/EpiStrainDynamics/issues/47)).
- Standardised error handling on
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
  ([\#49](https://github.com/acefa-hubs/EpiStrainDynamics/issues/49)).
- [`plot()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/plot.md)
  is no longer redefined as a generic
  ([\#48](https://github.com/acefa-hubs/EpiStrainDynamics/issues/48)).
- Increased the minimum R version to 4.1.0, matching actual base pipe
  usage, and added it to the R-CMD-check.yaml test matrix
  ([\#43](https://github.com/acefa-hubs/EpiStrainDynamics/issues/43)).
- Fixed all lintr-identified issues
  ([\#50](https://github.com/acefa-hubs/EpiStrainDynamics/issues/50)):
  styling (line length, indentation, spacing) via `styler`, dropped
  explicit [`return()`](https://rdrr.io/r/base/function.html) on final
  statements, converted [`sapply()`](https://rdrr.io/r/base/lapply.html)
  to [`vapply()`](https://rdrr.io/r/base/lapply.html) for type-safe
  extraction, and removed genuine dead code flagged by unused-variable
  checks.

### Documentation

- Replaced `\dontrun{}` in examples with `@examplesIf interactive()` for
  slow examples and guards such as
  `@examplesIf rlang::is_installed("xts")` where a package is used only
  in an example
  ([\#36](https://github.com/acefa-hubs/EpiStrainDynamics/issues/36)).
- Expanded documentation for
  [`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md),
  including its `gi_dist` argument
  ([\#44](https://github.com/acefa-hubs/EpiStrainDynamics/issues/44)).
- Expanded documentation for
  [`smoothing_structure()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/smoothing_structure.md)
  and
  [`diagnose_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/diagnose_model.md)
  ([\#53](https://github.com/acefa-hubs/EpiStrainDynamics/issues/53),
  [\#54](https://github.com/acefa-hubs/EpiStrainDynamics/issues/54)).
- Expanded documentation for
  [`growth_rate()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/growth_rate.md),
  [`incidence()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/incidence.md),
  and
  [`proportion()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/proportion.md)
  to match the level of detail in
  [`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md),
  including symbol breakdowns for their adjustment formulas and, for
  [`growth_rate()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/growth_rate.md),
  its relationship to
  [`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md)
  ([\#66](https://github.com/acefa-hubs/EpiStrainDynamics/issues/66)).
- Added roxygen `[]` cross-links for referenced functions throughout the
  documentation
  ([\#55](https://github.com/acefa-hubs/EpiStrainDynamics/issues/55)).
- Committed `precompile.R` for the vignette and included revision
  instructions in CONTRIBUTING.md
  ([\#46](https://github.com/acefa-hubs/EpiStrainDynamics/issues/46)).
- Clarified contribution and maintenance guidelines, added installation-
  from-source instructions, and added citations to the README
  ([\#51](https://github.com/acefa-hubs/EpiStrainDynamics/issues/51),
  [\#57](https://github.com/acefa-hubs/EpiStrainDynamics/issues/57),
  [\#58](https://github.com/acefa-hubs/EpiStrainDynamics/issues/58),
  [\#59](https://github.com/acefa-hubs/EpiStrainDynamics/issues/59),
  [\#60](https://github.com/acefa-hubs/EpiStrainDynamics/issues/60)).
- Corrected documentation of both `sarscov2` and `influenza` datasets
  ([\#56](https://github.com/acefa-hubs/EpiStrainDynamics/issues/56)).
- Added explanations of “pathogen noise” and the day-of-week effect
  ([\#61](https://github.com/acefa-hubs/EpiStrainDynamics/issues/61)), a
  link to the relevant timeseries class
  ([\#62](https://github.com/acefa-hubs/EpiStrainDynamics/issues/62)),
  and a COVID-19 example
  ([\#63](https://github.com/acefa-hubs/EpiStrainDynamics/issues/63)),
  to the vignette.
- Clarified the G1.6 srr standard statement to note that, although
  performance/scaling is evaluated in the vignettes, this does not
  constitute an explicit comparison against alternative implementations
  ([\#52](https://github.com/acefa-hubs/EpiStrainDynamics/issues/52)).

## EpiStrainDynamics 0.0.1 (2026-03-27)

Initial version submitted to rOpenSci for peer review
(ropensci/software-review#763).

### New features

- Random-walk and p-spline models for inferring temporal trends of
  single and multiple (including subtyped) pathogens, fitted via Stan.
- Model construction and fitting via
  [`construct_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/construct_model.md)
  and
  [`fit_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/fit_model.md),
  with
  [`smoothing_structure()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/smoothing_structure.md)
  and
  [`dispersion_structure()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/dispersion_structure.md)
  helpers.
- Post-processing generics
  [`incidence()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/incidence.md),
  [`growth_rate()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/growth_rate.md),
  [`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md),
  and
  [`proportion()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/proportion.md),
  each with a
  [`plot()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/plot.md)
  method.
- [`diagnose_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/diagnose_model.md)
  for MCMC convergence diagnostics.
- Bundled `influenza` and `sarscov2` example datasets.
- “Using EpiStrainDynamics” and “Algorithmic Scaling” vignettes.
