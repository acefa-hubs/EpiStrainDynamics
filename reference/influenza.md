# World Health Organisation Global Influenza Programme for Australia

Influenza-like illness data were retrieved from the World Health
Organization’s Global Influenza Programme. We collated weekly data for
Australia from the week starting January 2, 2012 to the week starting
December 25, 2023 inclusive. The data described: (1) the weekly number
of cases of influenza-like illness; and (2) the weekly number of
specimens positive for influenza by subtype. We grouped the influenza
specimens into: influenza A subtype not determined; influenza A H3N2;
influenza A H1N1 (influenza A H1N1, influenza A H1N1pdm09); influenza B
(influenza B Yamagata, influenza B Victoria, influenza B lineage not
determined); and other (unknown causes which we do not have data to
determine).

## Usage

``` r
influenza
```

## Format

### `influenza`

A data frame with 830 rows and 15 columns:

- ili:

  Integer, daily total number of cases of influenza-like illness

- week:

  Date, between 1 January 2012 to week starting 1 March 2020

- inf_A:

  Integer, daily number of cases of unsubtyped influenza A

- inf_B:

  Integer, daily number of cases of influenza B

- inf_H3N2:

  Integer, daily number of cases of influenza A subtype H3N2

- inf_H1N1:

  Integer, daily number of cases of influenza A subtype H1N1

- other:

  Integer, number of cases of unspecified influenza-like illness

## Source

<https://www.who.int/teams/global-influenza-programme/surveillance-and-monitoring/influenza-surveillance-outputs>
