#' United Kingdom Health Security Agency SARS-CoV-2 case data
#'
#' Daily SARS-CoV-2 case numbers by specimen date for the United Kingdom from
#' 2020 - 2022 retrieved from the UK Health Security Agency's data dashboard.
#' downloaded data describing the daily number of variants detected by
#' collection date. The data classified all sequences based on ‘major lineage
#' calls’. We considered 11 groupings of the variants based on their major
#' lineage calls or other, consisting of all lineages with a designation not
#' consistent with any of the major lineage calls.
#'
#' @format ## `sarscov2`
#' A data frame with 830 rows and 15 columns:
#' \describe{
#'   \item{date}{Date, between 23 September 2020 and 31 December 2022}
#'   \item{cases}{Numeric, daily number of cases}
#'   \item{B.1.177}{Integer, daily number of cases of B.1.177}
#'   \item{B.1.1.7}{Integer, daily number of cases of B.1.1.7 (Alpha variant)}
#'   \item{B.1.617.2}{Integer, daily number of cases of B.1.617.2 (Delta variant)}
#'   \item{BA.1}{Integer, daily number of cases of BA.1 (Omicron BA.1 variant)}
#'   \item{BA.2}{Integer, daily number of cases of BA.2 (Omicron BA.2 variant)}
#'   \item{BA.2.75}{Integer, daily number of cases of BA.4 (Omicron BA.4 variant)}
#'   \item{BA.4}{Integer, daily number of cases of BA.5 (Omicron BA.5 variant)}
#'   \item{BA.5}{Integer, daily number of cases of BA.2.75 (Omicron BA.2.75 variant)}
#'   \item{BQ.1}{Integer, daily number of cases of BQ.1 (Omicron BQ.1 variant)}
#'   \item{XBB}{Integer, daily number of cases of XBB (a recombinant of omicron
#'      sub-variants)}
#'   \item{Other}{Integer, daily number of cases of all lineages with a
#'      designation not consistent with any of the major lineage calls}
#'   \item{total}{}
#'   \item{t}{}
#' }
#' @source <https://ukhsa-dashboard.data.gov.uk/covid-19-archive-data-download;
#'   https://datadryad.org/dataset/doi:10.5061/dryad.hx3ffbgm2>
#' @srrstats {G5.1} Dataset `sarscov2` is exported and described in detail here.
#'   It is used in both example code and in testing.
#' @srrstats {G1.4} uses `Roxygen2` documentation
"sarscov2"

#' THIS NEEDS TO BE EDITED
#'
#' Influenza-like illness data were retrieved from the World Health
#' Organization’s Global Influenza Programme.16 We collated weekly data for
#' Australia from the week starting January 2, 2012 to the week starting
#' December 25, 2023 inclusive. The data described: (1) the weekly number
#' of cases of influenza-like illness; and (2) the weekly number of specimens
#' positive for influenza by subtype (and the number of negative tests). We
#' grouped the influenza specimens into: influenza A subtype not determined;
#' influenza A H3N2; influenza A H1N1 (influenza A H1N1, influenza A H1N1pdm09);
#' influenza B (influenza B Yamagata, influenza B Victoria, influenza B lineage
#' not determined); and negative tests. Negative tests likely reflect
#' individuals infected with other ILI-causing pathogens or other unknown
#' causes which we do not have data to determine.
#'
#' @format ## `influenza`
#' A data frame with 830 rows and 15 columns:
#' \describe{
#'   \item{ili}{Integer, daily total number of cases of influenza-like illness}
#'   \item{week}{Date, between 1 January 2012 to week starting 1 March 2020}
#'   \item{inf_A}{Integer, daily number of cases of unsubtyped influenza A}
#'   \item{inf_B}{Integer, daily number of cases of influenza B}
#'   \item{inf_all}{Integer, daily number of cases of }
#'   \item{inf_neg}{Integer, daily number of cases of }
#'   \item{inf_H3N2}{Integer, daily number of cases of influenza A subtype H3N2}
#'   \item{inf_H1N1}{Integer, daily number of cases of influenza A subtype H1N1}
#'   \item{num_spec}{Integer, number of specimens} ??
#' }
#' @source <https://www.who.int/teams/global-influenza-programme/surveillance-and-monitoring/influenza-surveillance-outputs>
#' @srrstats {G5.1} Dataset `influenza` is exported and described in detail
#'   here. It is used in both example code and in testing.
#' @srrstats {G1.4} uses `Roxygen2` documentation
"influenza"
