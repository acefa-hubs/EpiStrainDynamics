# Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data.

This repository recreates all analyses in the preprint:  Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data.

There are five analysis scripts in this directory that when ran will recreate all figures and supplementary figures from the preprint.

For each analysis script:
- The first (broad) section of code reads in and cleans publicly available data (found in the 'Data' subdirectory). 
- The second (broad) section of code then fits a series of stan models (this is the most time consuming component of the code) and then saves the fitted stan models as .RDS objects in the 'FitStanModels' subdirectory.
- The third (broad) section of code then reads in the fitted stan models (from the 'FitStanModels' subdirectory) and produces figures that are then saved to the 'Figures' subdirectory

## Influenza analysis
There are two scripts for producing the influenza analyses from the preprint.'influenza_analysis.R' produces all the main analysis and some of the supplementary analysis. 'influenza_sampling_rate_analysis.R' produces only some supplmentary analysis. 'influenza_analysis.R' needs to be run before 'influenza_sampling_rate_analysis.R' as the latter script relies on a stan model that is produced by the first script and saved in the 'FitStanModels' subdirectory.

## Dengue analyis
There are two scripts for producing the dengue analyses from the preprint.'dengue_analysis.R' produces all the main analysis and some of the supplementary analysis. 'dengue_loo_analysis.R' produces only some supplmentary analysis. 'dengue_analysis.R' needs to be run before 'dengue_loo_analysis.R' as the latter script relies on a stan model that is produced by the first script and saved in the 'FitStanModels' subdirectory.

## SARS-CoV-2 analysis
There is only one script that produces all SARS-CoV-2 analysis from the preprint: sarscov2_analysis.R
