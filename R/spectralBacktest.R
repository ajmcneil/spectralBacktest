#' spectralBacktest: Package for spectral backtesting of forecast distributions
#'
#' The spectralBacktest package implements spectral backtests of Z-test form as
#' defined in \href{https://arxiv.org/abs/1708.01489}{Gordy and McNeil (2018)}.  This package accommodates both
#' discrete and continuous kernels, and tests of both unconditional and conditional
#' coverage.
#'
#'
#' @section Monokernel, bikernel, multikernel and tlsfkernel objects:
#'
#' A monokernel is a list with the following named elements:
#' \itemize{
#'   \item{\code{name} is an optional name.}
#'   \item{\code{type} is string value 'mono'.}
#'   \item{\code{nu} is a closure taking parameters \code{support}, \code{param}, and \code{standardize},
#'     and returning a function mapping PIT values to transformed PIT values, i.e., \eqn{W=\nu(P)}.
#'     When \code{standardize} is \code{TRUE}, \eqn{W} will be standardized to mean zero and
#'     variance one under the null.}
#'   \item{\code{support} is the support of the kernel. For continuous kernels, \code{support} is
#'         a vector \eqn{(\alpha_1,\alpha_2)}.  For discrete kernels, \code{support} is the vector of
#'         points receiving positive weight.}
#'   \item{\code{param} is an optional list of parameters to the kernel,
#'          e.g., \eqn{\zeta} for the exponential kernel. For discrete kernels, \code{param} is the
#'          vector of weights on the points in \code{support}.}
#' }
#'
#' \preformatted{
#' Here is an example of a continuous monokernel:
#'      ZV <- list( name = 'Epanechnikov',
#'                  type = 'mono',
#'                  nu = nu_epanechnikov,
#'                  support = c(alpha1, alpha2),
#'                  param = NULL )
#' and an example of a discrete monokernel:
#'     ZU3 <- list( name = 'Discrete Uniform 3',
#'                  type = 'mono',
#'                  nu = nu_discrete,
#'                  support = c(alpha1, 0.99, alpha2),
#'                  param = c(1, 1, 1) )
#' The binomial score test is represented as:
#'     B99 <- list( name = 'Binomial score at 99pct',
#'                  type = 'mono',
#'                  nu = nu_discrete,
#'                  support = 0.99,
#'                  param = 1 )
#' }
#'
#' When writing your own monokernels (say, 'nu_custom'), it is a good idea to define a separate
#' function (say, 'mu_custom') to provide the first and second moments of the unstandardized kernel.
#' This can be useful when combining kernel functions into bikernel or multikernel objects.
#'
#' A bikernel is a list with the following named elements:
#' \itemize{
#'   \item{\code{name} is an optional name.}
#'   \item{\code{type} is string value 'bi'.}
#'   \item{\code{nu} is a list of closures, each as described in the monokernel case.}
#'   \item{\code{correlation} is a function taking parameters \code{support} and \code{param}, and returning
#'         the correlation between standardized transformed PIT values \eqn{W_1} and \eqn{W_2}.}
#'   \item{\code{support} is the common support of the kernels, as in the monokernel case.}
#'   \item{\code{param} is a list of optional parameter lists. Element \code{j} of list is passed to
#'          to the kernel in position \code{j} in \code{nu}. The full list is passed to the correlation
#'          function.}
#' }
#'
#' \preformatted{
#' Here is an example of a continuous bikernel:
#'     ZAE <- list( name = 'Arcsin/Epanechnikov',
#'                  type = 'bi',
#'                  nu = list(nu_arcsin, nu_epanechnikov),
#'                  correlation = rho_arcsin_epanechnikov,
#'                  support = c(alpha1, alpha2),
#'                  param = list(NULL, NULL) )
#' }
#'
#' For tractability, in the higher-order multikernel tests we assume that the component kernels are all in the
#' same family.  The simplifies the representation of a multikernel as a list with the following named elements:
#' \itemize{
#'   \item{\code{name} is a (string) name for the kernel.}
#'   \item{\code{type} is string value 'multi'.}
#'   \item{\code{nu} is a closure as in the monokernel case.}
#'   \item{\code{correlation} is a bivariate correlation function as in the bikernel case.}
#'   \item{\code{support} is the common support of the kernels as in the monokernel case.}
#'   \item{\code{param} is a list of optional parameter lists. Element \code{j} of list is passed to
#'          to the kernel in position \code{j} in \code{nu}. Pairs of elements are passed to the correlation
#'          function.}
#' }
#'
#' The Pearson multinomial test is expressed as a multikernel.  The component kernels do not
#' share a common support, but rather each has its own "level," and therefore the levels
#' are passed as a \code{param} list rather than in the \code{support} field:
#'  \preformatted{
#'     Pearson3 <- list( name = 'Pearson',
#'                       type = 'multi',
#'                       nu = nu_pearson,
#'                       correlation = rho_pearson,
#'                       support = NULL,
#'                       param = list(0.985, 0.99, 0.995) )
#' }
#'
#' Finally, we represent the TLSF bikernels as a list with the following named elements:
#' \itemize{
#'   \item{\code{name} is a (string) name for the kernel.}
#'   \item{\code{type} is string value 'tlsf'.}
#'   \item{\code{nu} is a closure taking parameters \code{support} and \code{param},
#'     and returning a function mapping PIT values to transformed PIT values, i.e., \eqn{W=\nu(P)}.
#'     \eqn{W} is a list in which element \eqn{j} is the j-th kernel applied to PIT. It is
#'     assumed that the kernels de-mean the transformed PIT but do not scale to unit variance.}
#'   \item{\code{VCV} is a function taking parameters \code{support} and \code{param},
#'     and returning a variance-covariance for the kernels under the null hypothesis.}
#'   \item{\code{support} is the support of the kernel represented as a vector \eqn{(\alpha_1,\alpha_2)}.}
#'   \item{\code{param} is an optional list of parameters to the kernel. Generally not used.}
#' }
#'
#' @section Tests of unconditional coverage:
#' Calls to the unconditional tests take a generic form:
#'
#' \preformatted{
#'   p_value <- spectral_Ztest(kernel, PIT, twosided)
#' }
#'
#' The \code{type} field in \code{kernel} redirects to specific handlers for mono/bi/multi.  For the bi- and multi-spectral cases, the
#' \code{twosided} parameter is unnecessary and would be ignored if supplied.
#'
#' @section Martingale difference tests:
#'
#' The conditional tests require an additional list structure to specify the \eqn{h}
#' function which is to be applied to the lagged PIT values. We refer to this structure as a
#' CVT (conditioning variable transformation).  This structure takes slightly different
#' forms in the monospectral and bispectral cases.
#'
#' A monospectral CVT is a list with the following elements:
#' \itemize{
#'   \item{\code{name} is a (string) name for the CVT.}
#'   \item{\code{type} is string value 'mono'.}
#'   \item{\code{h_closure} is a closure taking parameter \code{param}, and returning
#'         a function \eqn{h} mapping \eqn{P} values to \eqn{h(P)}.  Without loss of
#'         generality, it is required
#'         that the resulting function be standardized such that \eqn{h(P)} has mean
#'         zero and variance one when \eqn{P} is distributed Uniform.}
#'   \item{\code{h_param} is a list of optional parameters to \code{h_closure}.}
#'   \item{\code{lags} (integer) is the number of lagged PIT to include in \eqn{h}.}
#' }
#'
#' Here is an example based on the function \eqn{h(P) = |2P-1|^\theta} for \eqn{\theta>0}
#' \preformatted{
#'     CVT <- list( name = 'Power 4',
#'                  type = 'mono',
#'                  h_closure = CVT_PtildePower,
#'                  h_param = 4,
#'                  lags = 4L )
#' }
#' The closure \code{CVT_PtildePower} is included in the package.
#'
#' A bispectral CVT is a list with the following elements:
#' \itemize{
#'   \item{\code{name} is a (string) name for the CVT.}
#'   \item{\code{type} is string value 'bi'.}
#'   \item{\code{h_closure} is a list of two closures, each as in the monokernel case.}
#'   \item{\code{h_param} is a list of optional parameter lists. Element \code{j} of list is passed to
#'          to \code{h_closure[[j]]}.}
#'   \item{\code{lags} (vector of integers) is the number of lagged PIT to include in each \eqn{h}.}
#'   \item{\code{correlation} is a scalar function mapping \code{h_param} to the correlation of
#'         \code{h_closure[[1]]} and \code{h_closure[[2]]} under the null.}
#' }
#'
#' Here is an example with 4 lags and 2 lags applied to \eqn{\nu_1} and \eqn{\nu_2}, respectively:
#' \preformatted{
#'     CVT2 <- list( name = 'Power 4/Power 0.5',
#'                   type = 'bi',
#'                   h_closure = list(CVT_PtildePower, CVT_PtildePower),
#'                   correlation = rho_PtildePower,
#'                   h_param = list(4,1/2),
#'                   lags = c(4L,2L) )
#' }
#'
#' In principle, we can extend the MD tests to the multispectral case, but leave this for
#' a future version of the package.
#'
#' The MD test is called as follows:
#' \preformatted{
#'   p_value <- spectral_MDtest(kernel, CVT, PIT, h_list=NULL, estimateH=TRUE, contingencyH=FALSE)
#' }
#' where
#' \describe{
#'   \item{kernel}{ is a mono- or bi-kernel object as used in the unconditional tests.}
#'   \item{CVT}{ is a mono- or bi-spectral CVT object.}
#'   \item{PIT}{ is a vector of PIT values.}
#'   \item{h_list}{ (optional) is a precalculated matrix (in the monospectral case) or list of matrices (in the bispectral case).}
#'   \item{estimateH}{ (boolean). If \code{TRUE}, then estimated \eqn{H} is
#'      preferred whenever non-singular.  If \code{FALSE}, use the theoretical \eqn{H},
#'      which is the identity matrix in the monospectral case.}
#'   \item{contingencyH}{ (boolean).  If \code{TRUE}, then use theoretical \eqn{H}
#'        when the preferred \eqn{H} is singular. If \code{FALSE}, then return
#'        \code{NA} when preferred \eqn{H} singular.}
#' }
#'
#' The parameter \code{h_list} needs elaboration.  If this input is not provided by the user, the matrix (or list of matrices)
#' \eqn{h} are computed simply by applying the output of \code{h_closure} to the PIT values.  In the context of simulation studies,
#' however, it is often much more efficient to precalculate the \eqn{h} matrix for a set of tests. Allowing for \code{h_list} as
#' a specified input also permits flexibility in the situation in which inputed PIT values are substituted for missing PIT values
#' in generating \eqn{h}.
#'
#' @docType package
#' @name spectralBacktest
NULL

# Helper function to form a list of transformed PIT given a kernel object.
#   In monokernel case, returns a vector W rather than a list of W.
apply_nu_to_PIT <- function(kernel, PIT, standardize=TRUE) {
  if (kernel$type=='mono')
    return(kernel$nu(kernel$support, kernel$param, standardize)(PIT))

  if (kernel$type=='tlsf')
    return(kernel$nu(kernel$support, kernel$param)(PIT))

  isscaling <- 'scaling' == names(kernel$param)
  hasscaling <- any(isscaling)
  if (hasscaling) {
    stopifnot(length(kernel$param)==2)
    param <- kernel$param[!isscaling][[1]]
    scaling <- kernel$param[isscaling]$scaling
  } else {
    param <- kernel$param
  }

  # bikernel case
  if (kernel$type=='bi') {
    nuPIT <- function(k) {
      if (hasscaling) {
        prmsca <- list(param[[k]], scaling)
        names(prmsca) <- names(kernel$param)
        kernel$nu[[k]](kernel$support, prmsca, standardize)(PIT)
      } else {
        kernel$nu[[k]](kernel$support, param[[k]], standardize)(PIT)
      }
    }
    return(list(nuPIT(1), nuPIT(2)))
  }

  # multikernel case
  if (kernel$type=='multi') {
    nuPIT <- function(prm) {
      if (hasscaling) {
        prmsca <- list(prm, scaling)
        names(prmsca) <- names(kernel$param)
        kernel$nu(kernel$support, prmsca, standardize)(PIT)
      } else {
        kernel$nu(kernel$support, prm, standardize)(PIT)
      }
    }
    return(lapply(param, nuPIT))
  }
}

# Helper function to check for degeneracy in W.
#   In monokernel case, returns a vector W rather than a list of W.
is.degenerateW <- function(W) {
  if (is.list(W)) {
    purrr::map_lgl(W, function(ww) all(ww[1]==ww))
  } else {
    all(W[1]==W)
  }
}

# Helper function for exceedance.
#   Allows us to switch easily between GT and GTE conditions
exceedancePIT <- function(P,u)
  as.numeric(P>=u)

#' Theoretical correlation matrix for standardized multikernel
#'
#' @param multikernel (list) multikernel object (see package help)
#' @return \eqn{Sigma} correlation matrix of class Matrix
multispectral_correlation <- function(multikernel) {
  isscaling <- 'scaling' == names(multikernel$param)
  hasscaling <- any(isscaling)
  if (hasscaling) {
    stopifnot(length(multikernel$param)==2)
    param <- multikernel$param[!isscaling][[1]]
    scaling <- multikernel$param[isscaling]$scaling
  } else {
    param <- multikernel$param
  }
  J <- length(param)
  Sigma <- Matrix::Diagonal(J)
  rho <- function(jpair) {
    if (hasscaling) {
      multikernel$correlation(multikernel$support,
                              param[jpair],
                              scaling=scaling)
    } else {
      multikernel$correlation(multikernel$support, multikernel$param[jpair])
    }
  }

  for (j1 in (1:(J-1))) {
    for (j2 in ((j1+1):J)) {
      Sigma[j2,j1] <- Sigma[j1,j2] <- rho(c(j1,j2))
    }
  }
  return(Sigma)
}


#' Cross-moment to correlation.
#'
#' Returns a correlation given a cross-moment and two vectors of uncentered moments \eqn{(\mu_1,\mu_2)} for
#' each kernel.
#'
#' @param mux (double) cross-moment \eqn{E[\nu_A(P) \nu_B(P)]} under null
#' @param muA first two uncentered moments of \eqn{\nu_A(P)}
#' @param muB first two uncentered moments of \eqn{\nu_B(P)}
#' @param rho linear correlation of \eqn{\nu_A(P)} and \eqn{\nu_B(P)}
#' @export
cross_moment_to_correlation <- function(mux, muA, muB) {
  sigma <- sqrt((muA[2]-muA[1]^2)*(muB[2]-muB[1]^2))
  (mux - muA[1]*muB[1])/sigma |> unname()
}


#' Spectral Z-test for mono-kernel
#'
#' @param monokernel (list) mono-kernel object (see package help)
#' @param PIT Vector of probability integral transform values
#' @param twosided (boolean) run two-sided test if TRUE, otherwise one-sided
#' @return \eqn{p}-value of Z-test
#' @export
monospectral_Ztest <- function(monokernel, PIT, twosided=TRUE) {
  W <- PIT[!is.na(PIT)] |> as.numeric() |>
       monokernel$nu(monokernel$support, monokernel$param, standardize=TRUE)()
  if (is.degenerateW(W))
    warning(glue::glue('Kernel {monokernel$name}: transformed PIT has degenerate distribution.\n'))
  Z <- sqrt(length(W))*mean(W)
  if (twosided)
    output <- 2*(1-pnorm(abs(Z)))
  else
    output <- 1-pnorm(Z)
  return(output)
}

#' Spectral Z-test for bikernel
#'
#' @param bikernel (list) bikernel object (see package help)
#' @param PIT Vector of probability integral transform values
#' @return \eqn{p}-value of chisq-test
#' @export
bispectral_Ztest <- function(bikernel, PIT) {
  P <- as.numeric(PIT[!is.na(PIT)])
  # W_1 <- bikernel$nu[[1]](bikernel$support, bikernel$param[[1]], standardize=TRUE)(P)
  # W_2 <- bikernel$nu[[2]](bikernel$support, bikernel$param[[2]], standardize=TRUE)(P)
  # Wbar <- c(mean(W_1),mean(W_2))
  W_list <- apply_nu_to_PIT(bikernel, P, standardize = TRUE)
  if (any(is.degenerateW(W_list)))
    warning(glue::glue('Kernel {bikernel$name}: transformed PIT has degenerate distribution.\n'))
  Wbar <- purrr::map_dbl(W_list, mean)
#  rho <- bikernel$correlation(bikernel$support, bikernel$param)
#  Sigma <- matrix(c(1, rho, rho, 1), nrow=2)
  Sigma <- multispectral_correlation(bikernel)
  stat <-  length(P)*mahalanobis(Wbar, center=0, cov=Sigma)
  return(1-pchisq(stat,df=2))
}

#' Spectral Z-test for multikernel
#'
#' @param multikernel (list) multikernel object (see package help)
#' @param PIT Vector of probability integral transform values
#' @return \eqn{p}-value of chisq-test
#' @export
multispectral_Ztest <- function(multikernel, PIT) {
  P <- as.numeric(PIT[!is.na(PIT)])
  W_list <- apply_nu_to_PIT(multikernel, P, standardize = TRUE)
  if (any(is.degenerateW(W_list)))
    warning(glue::glue('Kernel {multikernel$name}: transformed PIT has degenerate distribution.\n'))
  Wbar <- purrr::map_dbl(W_list, mean)
  J <- length(Wbar)
  Sigma <- multispectral_correlation(multikernel)
  stat <-  length(P)*mahalanobis(Wbar, center=0, cov=Sigma)
  return(1-pchisq(stat,df=J))
}

#' Spectral Z-test for tlsfkernel
#'
#' @param tlsfkernel (list) TLSF kernel object (see package help)
#' @param PIT Vector of probability integral transform values
#' @return \eqn{p}-value of chisq-test
#' @export
tlsfspectral_Ztest <- function(tlsfkernel, PIT) {
  P <- na.omit(PIT)
  W_list <- apply_nu_to_PIT(tlsfkernel, P)
  if (any(is.degenerateW(W_list)))
    warning(glue::glue('Kernel {multikernel$name}: transformed PIT has degenerate distribution.\n'))
  Wbar <- purrr::map_dbl(W_list, mean)
  Sigma <- tlsfkernel$VCV(tlsfkernel$support, tlsfkernel$param)
  stat <-  length(P)*mahalanobis(Wbar, center=0, cov=Sigma)
  return(1-pchisq(stat,df=length(Wbar)))
}

#' Spectral Z-test for any kernel type
#'
#' @param kernel (list) kernel object (see package help)
#' @param PIT Vector of probability integral transform values
#' @param twosided (boolean) if monokernel, run two-sided test if TRUE, otherwise one-sided
#' @return \eqn{p}-value of Z-test
#' @export
spectral_Ztest <- function(kernel, PIT, twosided=TRUE) {
  switch(kernel$type,
         mono = monospectral_Ztest(kernel, PIT, twosided),
         bi = bispectral_Ztest(kernel, PIT),
         multi = multispectral_Ztest(kernel, PIT),
         tlsf = tlsfspectral_Ztest(kernel,PIT))
}
