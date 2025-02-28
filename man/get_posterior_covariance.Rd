% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Posterior_bands.R
\name{get_posterior_covariance}
\alias{get_posterior_covariance}
\title{Posterior inference for covariance surface}
\usage{
get_posterior_covariance(mcmc_results, zi, alpha_level = 0.05)
}
\arguments{
\item{mcmc_results}{mcmc object}

\item{zi}{Covariate vector of interest}

\item{alpha_level}{Type I error rate}
}
\value{
A cube
}
\description{
Posterior inference for covariance surface
}
