# Changelog

## EpiStrainDynamics (development version)

- Modified the plot code to explicitly call all ggplot2 functions and
  not import the whole package.
- Removed uses of in the examples. Slow examples use
  [@examplesIf](https://github.com/examplesIf) interactive() and other
  cases use alternatives, eg.
  [@examplesIf](https://github.com/examplesIf)
  rlang::is_installed(“xts”) where xts is used only in the example but
  not imported by the package.
- Standardised use of quotation marks style.
- Placed the arguments that do not have a default value before the ones
  that do.
- Added Config/testthat/parallel: true to DESCRIPTION file.
- Added additional tests for mcmc validations.

## EpiStrainDynamics 0.0.1.0000 (2026-03-27)

#### New features

- submitted to ROpenSci
