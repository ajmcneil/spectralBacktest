% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernelCollection.R
\name{truncatePIT}
\alias{truncatePIT}
\title{Truncate PIT to a kernel window and (optionally) rescale to unit interval.}
\usage{
truncatePIT(PIT, support = c(0, 1), rescale = FALSE)
}
\arguments{
\item{PIT}{Vector of probability integral transform values}

\item{support}{lower and upper bound of kernel window}
}
\value{
PITstar is vector of truncated PIT.
}
\description{
Calculate vector \eqn{P_t^* = \alpha_1 \vee (P_t \wedge \alpha_2)}.
If rescale==TRUE, then rescale to unit interval.
}
