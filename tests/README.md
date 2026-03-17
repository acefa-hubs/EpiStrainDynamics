# EpiStrainDynamics — test suite

## Running tests

Standard tests run with the usual `devtools::test()` or `R CMD check`. No
special configuration is needed. All tests, including regression tests against
pre-fitted models, run by default.

## How fixtures work

Rather than fitting Stan models during the test suite (which is slow), all
fitted model fixtures are pre-generated locally and stored in a GitHub Release.
They are downloaded automatically at the start of each test session via
[piggyback](https://docs.ropensci.org/piggyback/).

Two sets of fixtures are used:

- **Regular fixtures** (`helper-fixtures.R`): small data subsets, minimal MCMC
  settings (500 iterations, 1 chain). Used by `test-fit-model.R`,
  `test-diagnose-model.R`, `test-metrics.R`, and `test-plot.R`.
- **Extended fixtures** (`helper-extended-fixtures.R`): full package datasets,
  thorough MCMC settings (2000 iterations, 4 chains). Used by
  `test-extended-regression.R`.

Both helpers are sourced automatically by testthat before any tests run.
Fixtures live in memory for the duration of the test session and are discarded
afterwards.

## Regenerating fixtures

Fixtures only need to be regenerated when model structure changes in a way that
affects the fitted object format. To regenerate:
```r
devtools::load_all()
source("inst/scripts/generate-test-fixtures.R")
```

This fits all regular and extended fixtures locally, saves them to
`inst/testfixtures/`, and uploads them to the GitHub Release. Bump the tag
(e.g. `test-fixtures-v2`) and update the tag name in both helper files and the
generate script when doing so.

## Adding new tests

### Standard tests
Add test cases to the appropriate `test-*.R` file. If the test needs a new
fitted model, add it to `inst/scripts/generate-test-fixtures.R`, regenerate,
and load it in `helper-fixtures.R`.

### Extended regression tests
Add test cases to `test-extended-regression.R` (or a new `test-extended-*.R`
file). If the tests need a new fitted model, add it to the extended fixtures
section of `inst/scripts/generate-test-fixtures.R`, regenerate, and load it in
`helper-extended-fixtures.R`.

## Requirements

| Requirement | Detail |
|---|---|
| **Platform** | No platform-specific requirements. Tests run on Linux, macOS, and Windows. |
| **Memory** | Downloading and loading all fixtures requires roughly 200 MB of free RAM. |
| **Expected runtime** | Fast — no model fitting occurs during the test session. |
| **Artefacts** | No files are written to disk during tests. Fixtures are downloaded to a temporary directory and discarded afterwards. |
| **Dependencies** | Stan (via CmdStanR) must be installed to use fitted model objects. Run `cmdstanr::check_cmdstan()` to verify. |
