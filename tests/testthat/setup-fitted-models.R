# tests/testthat/helper-fixtures.R
# Downloads pre-fitted model fixtures from GitHub Release via piggyback.
# Fitting locally is skipped entirely on CI.

if (exists("sarscov2") && exists("influenza")) {

  fixture_dir <- file.path(tempdir(), "testfixtures")
  dir.create(fixture_dir, showWarnings = FALSE)

  message("Downloading pre-fitted model fixtures via piggyback...")

  piggyback::pb_download(
    repo    = "your-org/EpiStrainDynamics",
    tag     = "test-fixtures-v1",
    dest    = fixture_dir,
    overwrite = TRUE
  )

  # Load into global test environment
  fit_rw_single     <<- readRDS(file.path(fixture_dir, "fit_rw_single.rds"))
  fit_ps_single     <<- readRDS(file.path(fixture_dir, "fit_ps_single.rds"))
  fit_rw_multi      <<- readRDS(file.path(fixture_dir, "fit_rw_multi.rds"))
  fit_ps_multi      <<- readRDS(file.path(fixture_dir, "fit_ps_multi.rds"))
  fit_rw_subtyped   <<- readRDS(file.path(fixture_dir, "fit_rw_subtyped.rds"))
  fit_ps_subtyped   <<- readRDS(file.path(fixture_dir, "fit_ps_subtyped.rds"))
  fit_rw_single_dow <<- readRDS(file.path(fixture_dir, "fit_rw_single_dow.rds"))
  inc_single        <<- readRDS(file.path(fixture_dir, "inc_single.rds"))
  inc_multi         <<- readRDS(file.path(fixture_dir, "inc_multi.rds"))
  gr_single         <<- readRDS(file.path(fixture_dir, "gr_single.rds"))
  gr_multi          <<- readRDS(file.path(fixture_dir, "gr_multi.rds"))
  rt_single         <<- readRDS(file.path(fixture_dir, "rt_single.rds"))
  rt_multi          <<- readRDS(file.path(fixture_dir, "rt_multi.rds"))
  prop              <<- readRDS(file.path(fixture_dir, "prop.rds"))

  gi_simple <<- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))

  message("Fixtures ready.")

} else {
  message("Package data not available - skipping fixture setup")
}
