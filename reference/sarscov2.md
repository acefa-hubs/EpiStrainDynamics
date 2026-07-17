# United Kingdom Health Security Agency SARS-CoV-2 case data

Daily SARS-CoV-2 case numbers by specimen date for the United Kingdom
from 2020 - 2022 retrieved from the UK Health Security Agency's data
dashboard. downloaded data describing the daily number of variants
detected by collection date. The data classified all sequences based on
‘major lineage calls’. We considered groupings of the variants based on
their major lineage calls or other, consisting of all lineages with a
designation not consistent with any of the major lineage calls.

## Usage

``` r
sarscov2
```

## Format

### `sarscov2`

A data frame with 830 rows and 6 columns:

- date:

  Date, between 23 September 2020 and 31 December 2022

- cases:

  Numeric, daily number of cases

- alpha:

  Integer, daily number of cases of B.1.1.7 (Alpha variant)

- delta:

  Integer, daily number of cases of B.1.617.2 (Delta variant)

- omicron:

  Integer, daily number of cases of BA.1, BA.2, BA.2.75, BA.4, BA.5,
  BQ.1 (Omicron variants)

- other:

  Integer, daily number of cases of B.1.177, XBB (a recombinant of
  omicron sub-variants), and all other lineages with a designation not
  consistent with any of the major lineage calls

## Source

\<https://ukhsa-dashboard.data.gov.uk/covid-19-archive-data-download;
https://datadryad.org/dataset/doi:10.5061/dryad.hx3ffbgm2\>
