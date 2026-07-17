# Changelog

## EpiStrainDynamics (development version)

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
- Standardised error handling on
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
  ([\#49](https://github.com/acefa-hubs/EpiStrainDynamics/issues/49)).

### Documentation

- Replaced `\dontrun{}` in examples with `@examplesIf interactive()` for
  slow examples and guards such as
  `@examplesIf rlang::is_installed("xts")` where a package is used only
  in an example
  ([\#36](https://github.com/acefa-hubs/EpiStrainDynamics/issues/36)).
- Corrected documentation of both `sarscov2` and `influenza` datasets
  ([\#56](https://github.com/acefa-hubs/EpiStrainDynamics/issues/56)).
- Clarify contribution and maintenance guidelines, installation from
  source instructions, and add citations to readme
  ([\#51](https://github.com/acefa-hubs/EpiStrainDynamics/issues/51),
  [\#57](https://github.com/acefa-hubs/EpiStrainDynamics/issues/57),
  [\#58](https://github.com/acefa-hubs/EpiStrainDynamics/issues/58),
  [\#59](https://github.com/acefa-hubs/EpiStrainDynamics/issues/59),
  [\#60](https://github.com/acefa-hubs/EpiStrainDynamics/issues/60)).
- Expanded documentation for
  [`smoothing_structure()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/smoothing_structure.md)
  and
  [`diagnose_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/diagnose_model.md),
  add cross-linked references where applicable
  ([\#53](https://github.com/acefa-hubs/EpiStrainDynamics/issues/53),
  [\#54](https://github.com/acefa-hubs/EpiStrainDynamics/issues/54),
  [\#55](https://github.com/acefa-hubs/EpiStrainDynamics/issues/55)).
- Increased R version dependency in line with base pipe usage
  ([\#43](https://github.com/acefa-hubs/EpiStrainDynamics/issues/43)).

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
