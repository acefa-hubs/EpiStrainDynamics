# Pre-compiled vignettes that are slow to run
# Must manually move image files to vignettes/ after knit
knitr::knit("vignettes/Using-EpiStrainDynamics.Rmd.orig",
            output = "vignettes/Using-EpiStrainDynamics.Rmd")
