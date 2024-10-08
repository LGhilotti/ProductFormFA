% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GibbsModels.R
\name{GibbsFA_eb}
\alias{GibbsFA_eb}
\title{Gibbs-type feature allocation models through Empirical Bayes (GibbsFA_eb):
function to estimate the parameters via EB (maximizing the EFPF).}
\usage{
GibbsFA_eb(
  feature_matrix,
  model,
  type,
  seed = 1234,
  eb_params = NULL,
  Nhat_MM = NULL,
  var_fct = NULL,
  var_GammaIBP = NULL,
  ...
)
}
\arguments{
\item{feature_matrix}{A \code{n x K}-dimensional binary matrix of features}

\item{model}{Model to fit. Available models are \code{PoissonBB} (BB with Poisson(lambda) mixture),
\code{NegBinBB} (BB with NB(n0, mu0) mixture), \code{GammaIBP} (IBP with Gamma(a, b) mixture)}

\item{type}{Only "EFPF"}

\item{seed}{seed for fixing randomness}
}
\value{
An object of class \code{GibbsFA, model_eb}
}
\description{
Gibbs-type feature allocation models through Empirical Bayes (GibbsFA_eb):
function to estimate the parameters via EB (maximizing the EFPF).
}
