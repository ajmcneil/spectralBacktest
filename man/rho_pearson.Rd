% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernelCollection.R
\name{rho_pearson}
\alias{rho_pearson}
\title{Pearson correlation}
\usage{
rho_pearson(support = NULL, param = list(0.99, 0.99))
}
\arguments{
\item{support}{lower and upper bound of kernel window}

\item{param}{list of optional parameters to kernels}
}
\description{
Returns a bivariate correlation for VaR exceedance pairs for use in the Pearson test.
}
\details{
Note that \code{support} is ignored, and \code{param} takes its usual place.
This is to allow the Pearson test to be shoehorned into the multikernel
testing framework.on the desired \code{support} with weights in the \code{param} vector.
}
