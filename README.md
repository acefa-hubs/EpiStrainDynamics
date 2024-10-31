# Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data.

This repository contains code corresponding to the analyses in the preprint: _Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data_. 
For the R package in development based on these analyses, see branch [`main`](https://github.com/acefa-hubs/EpiStrainDynamics). 

Estimating the temporal trends in infectious disease activity is crucial for monitoring disease spread and the impact of interventions. 
Surveillance indicators routinely collected to monitor these trends are often a composite of multiple pathogens. 
For example, ‘influenza-like illness’ — routinely monitored as a proxy for influenza infections — is a symptom definition that could be caused by a wide range of pathogens, including multiple subtypes of influenza, SARS-CoV-2, and RSV. 
Inferred trends from such composite time series may not reflect the trends of any one of the component pathogens, each of which can exhibit distinct dynamics. 
Although many surveillance systems routinely test a subset of individuals contributing to a surveillance indicator — providing information on the relative contribution of the component pathogens — trends may be obscured by time-varying testing rates or substantial noise in the observation process. 
Here we develop a general statistical framework for inferring temporal trends of multiple pathogens from routinely collected surveillance data. 
We demonstrate its application to three different surveillance systems covering multiple pathogens (influenza, SARS-CoV-2, dengue), locations (Australia, Singapore, USA, Taiwan, UK), scenarios (seasonal epidemics, non-seasonal epidemics, pandemic emergence), and temporal reporting resolutions (weekly, daily). 
This methodology is applicable to a wide range of pathogens and surveillance systems.

## Analysis scripts
To reproduce this analysis you can clone the `preprint` branch of this repository as follows, replacing `<destination_folder>` with where you want to save the clone: 

```
git clone https://github.com/acefa-hubs/EpiStrainDynamics --branch preprint --single-branch <destination_folder>
```

There are five analysis scripts in this directory that will recreate all figures and supplementary figures from the preprint.

For each analysis script:
- The first (broad) section of code reads in and cleans publicly available data (found in the `Data` subdirectory). 
- The second (broad) section of code then fits a series of stan models (this is the most time consuming component of the code) and then saves the fitted stan models as .RDS objects in the `FitStanModels` subdirectory
- The third (broad) section of code then reads in the fitted stan models (from the `FitStanModels` subdirectory) and produces figures that are then saved to the `Figures` subdirectory

## Influenza analysis
There are two scripts for producing the influenza analyses from the preprint. 
`R/influenza_analysis.R` produces all the main analysis and some of the supplementary analysis. 
`R/influenza_sampling_rate_analysis.R` produces only some supplmentary analysis. 
`R/influenza_analysis.R` needs to be run before `R/influenza_sampling_rate_analysis.R` as the latter script relies on a stan model that is produced by the first script and saved in the `FitStanModels` subdirectory.

## Dengue analyis
There are two scripts for producing the dengue analyses from the preprint.
'dengue_analysis.R' produces all the main analysis and some of the supplementary analysis. 
`R/dengue_loo_analysis.R` produces only some supplmentary analysis. 
`R/dengue_analysis.R` needs to be run before `R/dengue_loo_analysis.R` as the latter script relies on a stan model that is produced by the first script and saved in the `FitStanModels` subdirectory.

## SARS-CoV-2 analysis
There is only one script that produces all SARS-CoV-2 analysis from the preprint: `R/sarscov2_analysis.R`

## Additional functions
Additional functions that are called by the main analysis scripts are housed in `R/inf_sampler.R`, `R/SFig12-13_functions.R`, and `R/SFig14_functions.R`. 

## Contribution
This is a work in progress. 
If you see any mistakes in the package (branch `main`) or in the code associated with the analyses for the preprint (branch `preprint`), let us know by logging a bug on the [issues](https://github.com/acefa-hubs/EpiStrainDynamics/issues) page. 
