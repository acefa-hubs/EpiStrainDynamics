% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_model_type.R
\name{get_noise_structure}
\alias{get_noise_structure}
\title{Get noise structure from specified observation noise}
\usage{
get_noise_structure(observation_noise)
}
\arguments{
\item{observation_noise}{either 'observation_noise_only' (only includes
observation noise - the same between pathogens) or pathogen_specific_noise'
(includes noise in individual pathogens as well)}
}
\value{
numeric value for noise structure needed by stan models
}
\description{
Get noise structure from specified observation noise
}
