#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.2} a Life Cycle Statement is available in `CONTRIBUTING.md`
#' @noRd
NULL

#' NA_standards
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#' @srrstatsNA {G1.6} This software makes no performance claims.
#' @srrstatsNA {G2.4c} no conversion to character object
#' @srrstatsNA {G2.4d, G2.4e, G2.5} no use of factor type objects
#' @srrstatsNA {G2.6} vector input is confined to columns names, so this
#'   pre-processing appears unnecessary
#' @srrstatsNA {G2.9} no information is lost from type conversion to tsibble
#' @srrstatsNA {G2.11, G2.12} data.frame-like objects are not used
#' @srrstatsNA {G2.14b, G2.14c} These models do not currently accommodate data
#'   with missing entries. So missing data cannot be ignored or replaced.
#' @srrstatsNA {G3.1, G3.1a} The software does not rely on covariance calculations
#' @srrstatsNA {G4.0} There are no outputs written to local files
#' @srrstatsNA {G5.0} This software is applicable to specific data so NIST
#'   datasets would not be applicable, and datasets are instead provided for
#'   example code and testing purposes.
#' @srrstatsNA {G5.4a} This software is not proposing a new algorithm and
#'   instead uses functions from the `rstan` package.
#' @srrstatsNA {G5.4b} This is not a new implementation of an existing method
#' @srrstatsNA {G5.4c} No stored values are needed for correctness testing
#' @srrstatsNA {G5.7} This software relies on Stan which has it's own underlying
#'  algortithm performance tests.
#' @srrstatsNA {G5.10,G5.11,G5.11a,G5.12} I have not implemented any extended
#'   tests
#' @srrstatsNA {BS1.0} The term 'hyperparameter' is not used
#' @srrstatsNA {BS1.3a, BS2.8} Since all the sampling functions use `sampling`
#'   from `rstan`, the user can specify initial values
#' @srrstatsNA {BS1.3b} This software uses Stan that implements only one type
#'   of algorithm (HMC/NUTS)
#' @srrstatsNA {BS1.5} This software provides a helper function to diagnose
#'   convergence, reporting the maximum R-hat and minimum n_eff. It does not
#'   explicitly choose for the user if the model converged, and so comparison
#'   of alternative dianostics is not necessary
#' @srrstatsNA {BS2.1, BS2.1a} only one tabular input, so numbers of rows don't
#'   need to be compared between multiple inputs
#' @srrstatsNA {BS2.6} Stan takes care of this directly
#' @srrstatsNA {BS2.9} This software uses `rstan`'s `sampling`, which ensures
#'   different seeds are used (from rstan's sampling() documentation: "Even if
#'   multiple chains are used, only one seed is needed, with other chains having
#'   seeds derived from that of the first chain to avoid dependent samples")
#' @srrstatsNA {BS2.10} This software ensures that the same seed is not passed
#'   to multiple computational chains
#' @srrstatsNA {BS2.13, BS2.14} This software uses `rstan`'s `sampling`, and
#'   while it seems possible to suppress all messages/warnings/progress bars, it
#'   is unclear how to suppress just progress bars or just warnings.
#' @srrstatsNA {BS3.1, BS3.2} These are not predictive models where there may be
#'   risk of collinear response and predictor variables
#' @srrstatsNA {BS4.0, BS4.1} This software does not implement internal
#'  sampling algorithms
#' @srrstatsNA {BS4.2, BS7.2} The only way to validate the posterior estimates would
#'  be to compare against an analytical solution, which we don't have.
#' @srrstatsNA {BS4.4} This software uses `rstan`'s `sampling`, and there does
#'   not appear to yet be a mechanism of stopping the chain upon convergence
#' @srrstatsNA {BS4.6} Convergence checker is a separate diagnostic function,
#'   it is not resulting in a shortened run sequence that needs to be evaluated
#' @srrstatsNA {BS4.7} This software does not implement a convergence checker
#'   whose threshold determines when to stop collecting samples
#' @srrstatsNA {BS5.4} A single function to diagnose convergence is used
#' @srrstatsNA {BS6.0} The fitting function, `fit_model()`, returns a list of
#'   class `EpiStrainDynamics.fit` as well as a class for the family of pathogen
#'   structure used (`ps`, `rw`, `ps_single`, or `rw_single`). The first item
#'   in the list is the stan output of class `stanfit`, which can be printed
#'   using default `rstan` print methods. By default the list is output to the
#'   console. A user can further interrogate the fit itself, or calculate
#'   epidemiological metrics using the `metrics` family of functions.
#' @srrstatsNA {BS6.4} `rstan` built-in summary methods can be used on the
#'   fitted model object
#' @srrstatsNA {BS7.0, BS7.1} Recovery of parametric estimates of the prior
#'   distributions are not scientifically relevant to this type of analysis.
#'   The priors are mainly implemented to improve sampling efficiency of the
#'   No-U-Turns Sampler in Stan and they are not related to the outcomes a user
#'   would be interested in evaluating.
#' @srrstatsNA {G5.6, G5.6a, G5.6b} Parameter recovery tests have not been
#'   implemented because the methods implemented in these stan models have
#'   been well tested previously in the literature. See README for full
#'   list of citations.
#' @noRd
NULL
