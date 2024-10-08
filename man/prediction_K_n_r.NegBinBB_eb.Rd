% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/K_n_r.R
\name{prediction_K_n_r.NegBinBB_eb}
\alias{prediction_K_n_r.NegBinBB_eb}
\title{Model-based K_n_r for BB with NB(n0,mu0) mixture - EB version}
\usage{
\method{prediction_K_n_r}{NegBinBB_eb}(object, n)
}
\arguments{
\item{object}{An object of class \code{GibbsFA, NegBinBB}}

\item{n}{Size of the hypothetical observed sample}
}
\description{
Model-based K_n_r for BB with NB(n0,mu0) mixture - EB version
}
\details{
Return list of the (n0,mu0)s sequences of the NB of K_n_r, for different n,
for BB with NB(n0,mu0) mixture
}
