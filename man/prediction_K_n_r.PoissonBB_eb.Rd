% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/K_n_r.R
\name{prediction_K_n_r.PoissonBB_eb}
\alias{prediction_K_n_r.PoissonBB_eb}
\title{Model-based K_n_r for BB with Poisson(lambda) mixture - EB version}
\usage{
\method{prediction_K_n_r}{PoissonBB_eb}(object, n)
}
\arguments{
\item{object}{An object of class \code{GibbsFA, PoissonBB_eb}}

\item{n}{Size of the hypothetical observed sample}
}
\description{
Model-based K_n_r for BB with Poisson(lambda) mixture - EB version
}
\details{
Return list of lambdas sequence of the K_n_r, for different n,
for BB with Poisson(lambda) mixture
}
