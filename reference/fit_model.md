# Generic Method for fitting model

S3 generic for fitted models from constructed model object

## Usage

``` r
fit_model(
  constructed_model,
  n_chain = 4,
  n_iter = 2000,
  n_warmup = floor(n_iter/2),
  thin = 1,
  adapt_delta = 0.9,
  multi_cores = TRUE,
  verbose = TRUE,
  suppress_warnings = FALSE,
  seed = NULL,
  ...
)

# S3 method for class 'rw_subtyped'
fit_model(
  constructed_model,
  n_chain = 4,
  n_iter = 2000,
  n_warmup = floor(n_iter/2),
  thin = 1,
  adapt_delta = 0.9,
  multi_cores = TRUE,
  verbose = TRUE,
  suppress_warnings = FALSE,
  seed = NULL,
  ...
)

# S3 method for class 'ps_subtyped'
fit_model(
  constructed_model,
  n_chain = 4,
  n_iter = 2000,
  n_warmup = floor(n_iter/2),
  thin = 1,
  adapt_delta = 0.9,
  multi_cores = TRUE,
  verbose = TRUE,
  suppress_warnings = FALSE,
  seed = NULL,
  ...
)

# S3 method for class 'rw_multiple'
fit_model(
  constructed_model,
  n_chain = 4,
  n_iter = 2000,
  n_warmup = floor(n_iter/2),
  thin = 1,
  adapt_delta = 0.9,
  multi_cores = TRUE,
  verbose = TRUE,
  suppress_warnings = FALSE,
  seed = NULL,
  ...
)

# S3 method for class 'ps_multiple'
fit_model(
  constructed_model,
  n_chain = 4,
  n_iter = 2000,
  n_warmup = floor(n_iter/2),
  thin = 1,
  adapt_delta = 0.9,
  multi_cores = TRUE,
  verbose = TRUE,
  suppress_warnings = FALSE,
  seed = NULL,
  ...
)

# S3 method for class 'rw_single'
fit_model(
  constructed_model,
  n_chain = 4,
  n_iter = 2000,
  n_warmup = floor(n_iter/2),
  thin = 1,
  adapt_delta = 0.9,
  multi_cores = TRUE,
  verbose = TRUE,
  suppress_warnings = FALSE,
  seed = NULL,
  ...
)

# S3 method for class 'ps_single'
fit_model(
  constructed_model,
  n_chain = 4,
  n_iter = 2000,
  n_warmup = floor(n_iter/2),
  thin = 1,
  adapt_delta = 0.9,
  multi_cores = TRUE,
  verbose = TRUE,
  suppress_warnings = FALSE,
  seed = NULL,
  ...
)
```

## Arguments

- constructed_model:

  prepared model object of class `EpiStrainDynamics.model`

- n_chain:

  number of MCMC chains, defaults to 4

- n_iter:

  A positive integer specifying the number of iterations for each chain,
  default value is 2000

- n_warmup:

  A positive integer specifying the number of warmup iterations, default
  value is half the number of iterations

- thin:

  A positive integer specifying the period for saving samples, default
  value is 1.

- adapt_delta:

  Numeric value between 0 and 1 indicating target average acceptance
  probability used in
  [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
  Default value is 0.9.

- multi_cores:

  A logical value indicating whether to parallelize chains with multiple
  cores, default is TRUE and uses all available cores - 1.

- verbose:

  Logical value controlling the verbosity of output. When TRUE
  (default), shows all messages, warnings, errors, and progress
  indicators. When FALSE, suppresses messages and progress while
  retaining warnings and errors.

- suppress_warnings:

  Logical value indicating whether to suppress warnings from Stan.
  Default is FALSE. When TRUE, warnings are suppressed but errors are
  still raised.

- seed:

  A positive integer seed used for random number generation in MCMC.
  Default is NULL, which means the seed is generated from 1 to the
  maximum integer supported by R.

- ...:

  additional arguments to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html),
  such as `init`

## Value

fit model of class `EpiStrainDynamics.fit`, or if fitting fails, an
error is raised that can be caught and inspected.

## Examples

``` r
if (FALSE) { # \dontrun{
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date))

  fit <- fit_model(mod)

  # Suppress progress and messages but keep warnings/errors
  fit <- fit_model(mod, verbose = FALSE)

  # Suppress warnings too
  fit <- fit_model(mod, verbose = FALSE, suppress_warnings = TRUE)

  # Catch errors and inspect
  result <- tryCatch(
    fit_model(mod),
    error = function(e) e
  )
  if (inherits(result, "EpiStrainDynamics.fit.error")) {
    cat("Fitting failed:", result$message, "\n")
    # Can still access the model: result$constructed_model
  }
} # }
```
