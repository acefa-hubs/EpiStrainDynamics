# Pre-compiled vignettes that are slow to run

devtools::load_all()

knitr::opts_knit$set(base.dir = "vignettes")
knitr::knit("vignettes/Using-EpiStrainDynamics.Rmd.orig",
            output = "vignettes/Using-EpiStrainDynamics.Rmd")

knitr::opts_knit$set(base.dir = "vignettes/articles")
knitr::knit("vignettes/articles/parameter-recovery.Rmd.orig",
            output = "vignettes/articles/parameter-recovery.Rmd")

knitr::opts_knit$set(base.dir = "vignettes/articles")
knitr::knit("vignettes/articles/algorithmic-scaling.Rmd.orig",
            output = "vignettes/articles/algorithmic-scaling.Rmd")
