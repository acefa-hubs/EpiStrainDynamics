# Pre-compiled vignettes that are slow to run

knitr::opts_knit$set(base.dir = "vignettes")
knitr::knit("vignettes/Using-EpiStrainDynamics.Rmd.orig",
            output = "vignettes/Using-EpiStrainDynamics.Rmd")
