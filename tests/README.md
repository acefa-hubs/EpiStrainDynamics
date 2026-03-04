# EpiStrainDynamics — test suite

## Running tests

Standard tests run with the usual `devtools::test()` or `R CMD check`. No
special configuration is needed.

## Extended tests

A subset of tests fit Stan models and are therefore slow (several hours on a
typical machine). These tests are **skipped by default** and are gated behind
an environment variable so that they can be opted into deliberately.

### Enabling extended tests locally

Set the environment variable **before** starting R (or export it in your shell):

  ```bash
export EPISTRAINDYNAMICS_EXTENDED_TESTS=true
Rscript -e "devtools::test()"
```

Or set it within an R session before calling `devtools::test()`:

```r
Sys.setenv(EPISTRAINDYNAMICS_EXTENDED_TESTS = "true")
devtools::test() # OR devtools::test(filter = "extended")
```

### How the extended tests work

1. `tests/testthat/helper-extended-fixtures.R` is sourced automatically by
testthat before any tests run. When the environment variable is set it fits
all fixture models (random walk and p-spline, single and multiple pathogen)
and stores the results in a shared in-memory environment called
`extended_fixtures`.

2. `tests/testthat/test-extended-regression.R` contains the actual test cases.
Each test calls `skip_extended()` at the top; if the environment variable is
not set the entire file is skipped with a single diagnostic message.

3. The fixtures use package datasets so are fully reproducible. No external 
data downloads are required.

### Requirements and expectations

| Requirement | Detail |
  |---|---|
  | **Platform** | No platform-specific requirements. Tests run on Linux, macOS, and Windows. |
  | **Memory** | Fitting all fixtures requires roughly 200 MB of free RAM. |
  | **Expected runtime** | around 4 hours depending on the machine. |
  | **Artefacts** | No files are written to disk. All fixtures live in memory for the duration of the test session and are discarded afterwards. |
  | **Dependencies** | Stan (via CmdStanR) must be installed. Run `cmdstanr::check_cmdstan()` to verify. |

### Adding new extended tests

1. Add test cases to `test-extended-regression.R` (or a new
`test-extended-*.R` file).
2. Gate every `test_that` block with `skip_extended()` line.
3. If the tests need a new fitted model, add it to `helper-extended-fixtures.R`
inside the existing `if (extended_tests_enabled())` block and store it on
`extended_fixtures`.
