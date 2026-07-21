#' United Kingdom Health Security Agency SARS-CoV-2 case data
#'
#' Daily SARS-CoV-2 case numbers by specimen date for the United Kingdom from
#' 2020 - 2022 retrieved from the UK Health Security Agency's data dashboard.
#' downloaded data describing the daily number of variants detected by
#' collection date. The data classified all sequences based on ‘major lineage
#' calls’. We considered groupings of the variants based on their major
#' lineage calls or other, consisting of all lineages with a designation not
#' consistent with any of the major lineage calls.
#'
#' @format ## `sarscov2`
#' A data frame with 830 rows and 6 columns:
#' \describe{
#'   \item{date}{Date, between 23 September 2020 and 31 December 2022}
#'   \item{cases}{Numeric, daily number of cases}
#'   \item{alpha}{Integer, daily number of cases of B.1.1.7 (Alpha variant)}
#'   \item{delta}{Integer, daily number of cases of B.1.617.2 (Delta variant)}
#'   \item{omicron}{Integer, daily number of cases of BA.1, BA.2, BA.2.75,
#'      BA.4, BA.5, BQ.1 (Omicron variants)}
#'   \item{other}{Integer, daily number of cases of B.1.177, XBB (a recombinant
#'      of omicron sub-variants), and all other lineages with a designation
#'      not consistent with any of the major lineage calls}
#' }
#' @source <https://ukhsa-dashboard.data.gov.uk/covid-19-archive-data-download;
#'   https://datadryad.org/dataset/doi:10.5061/dryad.hx3ffbgm2>
#' @srrstats {G5.1} Dataset `sarscov2` is exported and described in detail here.
#'   It is used in both example code and in testing.
#' @srrstats {G1.4} uses `Roxygen2` documentation
"sarscov2"

#' World Health Organisation Global Influenza Programme for Australia
#'
#' Influenza-like illness data were retrieved from the World Health
#' Organization’s Global Influenza Programme. We collated weekly data for
#' Australia from the week starting January 2, 2012 to the week starting
#' December 25, 2023 inclusive. The data described: (1) the weekly number
#' of cases of influenza-like illness; and (2) the weekly number of specimens
#' positive for influenza by subtype. We grouped the influenza specimens into:
#' influenza A subtype not determined; influenza A H3N2; influenza A H1N1
#' (influenza A H1N1, influenza A H1N1pdm09); influenza B (influenza B Yamagata,
#' influenza B Victoria, influenza B lineage not determined); and other (unknown
#' causes which we do not have data to determine).
#'
#' @format ## `influenza`
#' A data frame with 426 rows and 7 columns:
#' \describe{
#'   \item{ili}{Integer, daily total number of cases of influenza-like illness}
#'   \item{week}{Date, between 1 January 2012 to week starting 1 March 2020}
#'   \item{inf_A}{Integer, daily number of cases of unsubtyped influenza A}
#'   \item{inf_B}{Integer, daily number of cases of influenza B}
#'   \item{inf_H3N2}{Integer, daily number of cases of influenza A subtype H3N2}
#'   \item{inf_H1N1}{Integer, daily number of cases of influenza A subtype H1N1}
#'   \item{other}{Integer, number of cases of unspecified influenza-like illness}
#' }
#' @source <https://www.who.int/teams/global-influenza-programme/surveillance-and-monitoring/influenza-surveillance-outputs>
#' @srrstats {G5.1} Dataset `influenza` is exported and described in detail
#'   here. It is used in both example code and in testing.
#' @srrstats {G1.4} uses `Roxygen2` documentation
"influenza"
