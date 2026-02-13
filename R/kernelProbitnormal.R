# These functions have all been deprecated in favor of the generalized TLSF approach.

# #' Probitnormal score correlation
# #'
# #' Returns the correlation for VaR exceedance pairs for use in the PNS test.
# #'
# #' @inheritParams .bikernel_correlation
# #' @export
# rho_probitnormal_old <- function(support=c(0,1),param=NULL) {
#   qq <- qnorm(support)
#   dd <- dnorm(qq)
#   mux <- dd[1]*(1+qq[1]^2) - dd[2]*(1+qq[2]^2) +
#     (qq[1]*dd[1]^2)/support[1] + (qq[2]*dd[2]^2)/(1-support[2])
#   rho <-cross_moment_to_correlation(mux,
#                                     mu_probitnormal_old(support, param=1),
#                                     mu_probitnormal_old(support, param=2))
#   return(rho)
# }

# Probitnormal Score case

# #' Probitnormal score kernel1 function
# #'
# #' Returns a function nu_1(P), where \eqn{\nu_1} is the centered PNS kernel1 on
# #' the desired \code{support}.
# #'
# #' @inheritParams .monokernel_function
# #' @export
# nu_probitnormal1 <- function(support=c(0,1), param=NULL, standardize=TRUE) {
#   nu1 <- function(PIT) {
#     alpha1 <- support[1]
#     alpha2 <- support[2]
#     nu1_low <- -dnorm(qnorm(alpha1))/alpha1
#     nu1_high <- dnorm(qnorm(alpha2))/(1-alpha2)
#     score1 <- function(p) {
#       Pregion <- cut(p, breaks=c(0, alpha1, alpha2, 1),
#                      labels=c('low','mid','high'), include.lowest = TRUE)
#       switch(as.character(Pregion),
#              low = nu1_low,
#              mid = qnorm(p),
#              high = nu1_high )
#     }
#     W <- sapply(PIT,score1)
#     if (standardize) {
#       mu <- mu_probitnormal1(support)
#       W <- as.numeric((W-mu[1])/sqrt(mu[2]-mu[1]^2))
#     }
#     return(W)
#   }
#   return(nu1)
# }
#
# #' Probitnormal score kernel2 function
# #'
# #' Returns a function nu_2(P), where \eqn{\nu_2} is the centered PNS kernel2 on
# #' the desired \code{support}.
# #'
# #' @inheritParams .monokernel_function
# #' @export
# nu_probitnormal2 <- function(support=c(0,1), param=NULL, standardize=TRUE) {
#   nu2 <- function(PIT) {
#     alpha1 <- support[1]
#     alpha2 <- support[2]
#     nu2_low <- -qnorm(alpha1)*dnorm(qnorm(alpha1))/alpha1
#     nu2_high <- qnorm(alpha2)*dnorm(qnorm(alpha2))/(1-alpha2)
#     score2 <- function(p) {
#       Pregion <- cut(p, breaks=c(0, alpha1, alpha2, 1),
#                      labels=c('low','mid','high'), include.lowest = TRUE)
#       switch(as.character(Pregion),
#              low = nu2_low,
#              mid = qnorm(p)^2-1,
#              high = nu2_high )
#     }
#     W <- sapply(PIT,score2)
#     if (standardize) {
#       mu <- mu_probitnormal2(support)
#       W <- as.numeric((W-mu[1])/sqrt(mu[2]-mu[1]^2))
#     }
#     return(W)
#   }
#   return(nu2)
# }

# #' Probitnormal score kernel function
# #'
# #' Returns a function nu(P), where \eqn{\nu} is the centered PNS j-kernel (\eqn{j=1,2}) on
# #' the desired \code{support}. Set \code{param} to kernel \eqn{j}.
# #'
# #' @inheritParams .monokernel_function
# #' @export
# nu_probitnormal_old <- function(support=c(0,1), param=0, standardize=TRUE) {
#   nu <- function(PIT) {
#     stopifnot(param %in% c(1,2))
#     alpha1 <- support[1]
#     alpha2 <- support[2]
#     if (param==1) {
#       nu_low <- -dnorm(qnorm(alpha1))/alpha1
#       nu_high <- dnorm(qnorm(alpha2))/(1-alpha2)
#       mid_fn <- qnorm
#     } else {
#       nu_low <- -qnorm(alpha1)*dnorm(qnorm(alpha1))/alpha1
#       nu_high <- qnorm(alpha2)*dnorm(qnorm(alpha2))/(1-alpha2)
#       mid_fn <- function(p) qnorm(p)^2-1
#     }
#     scorefn <- function(p) {
#       stopifnot(p>=0, p<=1)
#       Pregion <- cut(p, breaks=c(0, alpha1, alpha2, 1.01),
#                      labels=c('low','mid','high'), include.lowest = TRUE)
#       switch(as.character(Pregion),
#              low = nu_low,
#              mid = mid_fn(p),
#              high = nu_high )
#     }
#     W <- map_dbl(PIT,scorefn)
#     if (standardize) {
#       mu <- mu_probitnormal_old(support, param) |> unname()
#       W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
#     }
#     return(W)
#   }
#   return(nu)
# }

#' Probitnormal score kernel1 moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' PNS kernel1 on the desired \code{support}.
#'
#' Deprecated function--use \code{nu_probitnormal} instead.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_probitnormal1 <- function(support=c(0,1), param=NULL) {
  qq <- qnorm(support)
  dd <- dnorm(qq)
  mu2 <- (dd[1]^2)/support[1] + (dd[2]^2)/(1-support[2]) +
    dd[1]*qq[1] - dd[2]*qq[2] + diff(support)
  mu <- c(mu1=0, mu2=mu2)
  return(mu)
}

#' Probitnormal score kernel2 moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' PNS kernel2 on the desired \code{support}.
#'
#' Deprecated function--use \code{nu_probitnormal} instead.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_probitnormal2 <- function(support=c(0,1), param=NULL) {
  qq <- qnorm(support)
  dd <- dnorm(qq)
  mu2 <- dd[1]*qq[1] - dd[2]*qq[2] + dd[1]*qq[1]^3 - dd[2]*qq[2]^3 +
    ((dd[1]*qq[1])^2)/support[1] + ((dd[2]*qq[2])^2)/(1-support[2]) +
    2*diff(support)
  mu <- c(mu1=0, mu2=mu2)
  return(mu)
}

# #' Probitnormal score kernel moments
# #'
# #' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
# #' PNS kernel on the desired \code{support}.  Set \code{param} to control kernel \eqn{j}.
# #'
# #' @inheritParams .monokernel_moments
# #' @export
# mu_probitnormal_old <- function(support=c(0,1), param=0) {
#   stopifnot(param %in% c(1,2))
#   qq <- qnorm(support)
#   dd <- dnorm(qq)
#   if (param==1) {
#     mu2 <- (dd[1]^2)/support[1] + (dd[2]^2)/(1-support[2]) +
#       dd[1]*qq[1] - dd[2]*qq[2] + diff(support)
#   } else {
#     mu2 <- dd[1]*qq[1] - dd[2]*qq[2] + dd[1]*qq[1]^3 - dd[2]*qq[2]^3 +
#       ((dd[1]*qq[1])^2)/support[1] + ((dd[2]*qq[2])^2)/(1-support[2]) +
#       2*diff(support)
#   }
#   mu <- c(mu1=0, mu2=mu2)
#   return(mu)
# }

#' Probitnormal score correlation
#'
#' Returns the correlation for VaR exceedance pairs for use in the PNS test.
#'
#' Deprecated function--use \code{vcv_probitnormal} instead.
#'
#' @inheritParams .bikernel_correlation
#' @export
rho_probitnormal <- function(support=c(0,1),param=NULL) {
  qq <- qnorm(support)
  dd <- dnorm(qq)
  mux <- dd[1]*(1+qq[1]^2) - dd[2]*(1+qq[2]^2) +
    (qq[1]*dd[1]^2)/support[1] + (qq[2]*dd[2]^2)/(1-support[2])
  rho <-cross_moment_to_correlation(mux,
                                    mu_probitnormal1(support),
                                    mu_probitnormal2(support))
  return(rho)
}

# rho_probitnormal_old <- function(support=c(0,1),param=NULL) {
#   qq <- qnorm(support)
#   dd <- dnorm(qq)
#   mux <- dd[1]*(1+qq[1]^2) - dd[2]*(1+qq[2]^2) +
#     (qq[1]*dd[1]^2)/support[1] + (qq[2]*dd[2]^2)/(1-support[2])
#   rho <-cross_moment_to_correlation(mux,
#                                     mu_probitnormal_old(support, param=1),
#                                     mu_probitnormal_old(support, param=2))
#   return(rho)
# }
