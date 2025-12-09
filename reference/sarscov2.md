# United Kingdom Health Security Agency SARS-CoV-2 case data

Daily SARS-CoV-2 case numbers by specimen date for the United Kingdom
from 2020 - 2022 retrieved from the UK Health Security Agency's data
dashboard. downloaded data describing the daily number of variants
detected by collection date. The data classified all sequences based on
‘major lineage calls’. We considered 11 groupings of the variants based
on their major lineage calls or other, consisting of all lineages with a
designation not consistent with any of the major lineage calls.

## Usage

``` r
sarscov2
```

## Format

### `sarscov2`

A data frame with 830 rows and 15 columns:

- date:

  Date, between 23 September 2020 and 31 December 2022

- cases:

  Numeric, daily number of cases

- B.1.177:

  Integer, daily number of cases of B.1.177

- B.1.1.7:

  Integer, daily number of cases of B.1.1.7 (Alpha variant)

- B.1.617.2:

  Integer, daily number of cases of B.1.617.2 (Delta variant)

- BA.1:

  Integer, daily number of cases of BA.1 (Omicron BA.1 variant)

- BA.2:

  Integer, daily number of cases of BA.2 (Omicron BA.2 variant)

- BA.2.75:

  Integer, daily number of cases of BA.4 (Omicron BA.4 variant)

- BA.4:

  Integer, daily number of cases of BA.5 (Omicron BA.5 variant)

- BA.5:

  Integer, daily number of cases of BA.2.75 (Omicron BA.2.75 variant)

- BQ.1:

  Integer, daily number of cases of BQ.1 (Omicron BQ.1 variant)

- XBB:

  Integer, daily number of cases of XBB (a recombinant of omicron
  sub-variants)

- Other:

  Integer, daily number of cases of all lineages with a designation not
  consistent with any of the major lineage calls

- total:

- t:

## Source

\<https://ukhsa-dashboard.data.gov.uk/covid-19-archive-data-download;
https://datadryad.org/dataset/doi:10.5061/dryad.hx3ffbgm2\>
