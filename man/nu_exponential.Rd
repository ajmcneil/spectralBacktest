% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernelCollection.R
\name{nu_exponential}
\alias{nu_exponential}
\title{Exponential kernel function}
\usage{
nu_exponential(support = c(0, 1), param = 0, standardize = TRUE)
}
\arguments{
\item{support}{lower and upper bound of kernel window}

\item{param}{optional parameters to kernel}

\item{standardize}{(boolean) if TRUE, then standardize \eqn{W} to mean 0, variance 1 under the null}
}
\description{
Returns a function nu(P), where nu is the exponential kernel on
the desired \code{support} with parameter \eqn{\zeta=}\code{param}.
}
