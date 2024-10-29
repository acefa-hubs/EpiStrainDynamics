### Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data.

This repository recreates all analyses in the preprint:  Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data.

There are five analysis scripts in this directory that when ran will recreate all figures and supplementary figures from the preprint.

For each analysis script:
- The first (broad) section of code reads in and cleans publicly available data (found in the 'Data' subdirectory). 
- The second (broad) section of code then fits a series of stan models (this is the most time consuming component of the code) and then saves the fitted stan models as .RDS objects in the 'FitStanModels' subdirectory.
- The third (broad) section of code then reads in the fitted stan models (from the 'FitStanModels' subdirectory) and produces figures that are then saved to the 'Figures' subdirectory
