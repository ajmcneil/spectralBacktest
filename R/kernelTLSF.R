#' Generic TLSF kernel function for inheritance of parameter documentation
#'
#' @param support lower and upper bound of kernel window
#' @param param optional parameters to kernel
.tlsfkernel_function <- function(support, param)
  NULL

#' VCV for TLSF
#'
#' Most of the heavy-lifting in computing \eqn{\Sigma} is here.
#' This function also checks that \eqn{\alpha_1} is above its admissible minimum.
#'
#' @details Current version of function assumes that the following terms go to
#' zero as \eqn{\alpha_2\rightarrow 1}:
#' * \code{alpha2+A0(alpha2)*Rinverse(alpha2)}
#' * \code{(1-alpha2)*CC2(alpha2,Cbar1)}
#' * \code{(1-alpha2)*(Cbar1(alpha2) + Rinverse(alpha2)*CC2(alpha2,Cbar1))}
#' * \code{(1-alpha2)*Rinverse(alpha2)*(2*Cbar1(alpha2) +
#'         Rinverse(alpha2)*CC2(alpha2,Cbar1))}
#' @md
#'
#' @param support lower and upper bound of kernel window
#' @param Rinverse quantile function of TLSF \eqn{R}
#' @param rho density function \eqn{\rho}
#' @param lambdarho function \eqn{\lambda_\rho}
#' @param B0,B1,B2 integrals expressed as functions of alpha
#' @return \eqn{I(\theta)} information matrix
#' @export
vcv_tlsf <- function(support, Rinverse, rho, lambdarho, B0, B1, B2) {
  alpha1 <- support[1]
  alpha2 <- support[2]
  C1 <- function(alf) rho(Rinverse(alf))/alf
  Cbar1 <- function(alf) -rho(Rinverse(alf))/(1-alf)
  gamma21<- Rinverse(alpha1)*(C1(alpha1)+lambdarho(Rinverse(alpha1)))-1
  stopifnot("Boundary alpha1 is out of range"=gamma21>=0,
            "Boundary alpha2 is out of range"=alpha2<=1)
  # CC2 accommodates C2 and Cbar2
  CC2 <- function(alf,CC) -CC(alf)*(lambdarho(Rinverse(alf))+CC(alf))
  A0 <- function(alf) if (abs(1-2*alf)==1) { 0 } else {-rho(Rinverse(alf))}
  A1 <- function(alf) if (abs(1-2*alf)==1) { alf } else {
           alf+A0(alf)*Rinverse(alf) }
  infomat <- matrix(0,nrow=2,ncol=2) # initialize
  if (alpha2<1) {
    tailpart11 <- (1-alpha2)*CC2(alpha2,Cbar1)
    tailpart12 <- (1-alpha2)*(Cbar1(alpha2)+Rinverse(alpha2)*CC2(alpha2,Cbar1))
    tailpart22 <- (1-alpha2)*Rinverse(alpha2)*(2*Cbar1(alpha2) +
                                     Rinverse(alpha2)*CC2(alpha2,Cbar1))
  } else {
    tailpart11 <- tailpart12 <- tailpart22 <- 0
  }
  infomat[1,1] <- B0(alpha2)-B0(alpha1) - alpha1*CC2(alpha1,C1) -
    tailpart11
  infomat[1,2] <- A0(alpha2)-A0(alpha1) + B1(alpha2)-B1(alpha1) -
    alpha1*(C1(alpha1) + Rinverse(alpha1)*CC2(alpha1,C1)) -
    tailpart12
  infomat[2,1] <- infomat[1,2]
  infomat[2,2] <- 2*(A1(alpha2)-A1(alpha1)) + B2(alpha2)-B2(alpha1) -
    (alpha2-alpha1) - alpha1*Rinverse(alpha1)*(2*C1(alpha1) +
                         Rinverse(alpha1)*CC2(alpha1,C1)) - tailpart22
  return(infomat)
}

#' Score kernel function
#'
#' Returns a function nu(P), where \eqn{\nu} itself returns a list with the
#' centered PNS kernels on the desired \code{support}.
#'
#' @param support lower and upper bound of kernel window
#' @param Rinverse quantile function of TLSF \eqn{R}
#' @param rho density function \eqn{\rho}
#' @param lambdarho function \eqn{\lambda_\rho}
#' @return \eqn{\nu(P)} kernel function
#' @export
tlsf_scorefunction <- function(support=c(0,1), Rinverse, rho, lambdarho) {
  nu <- function(PIT) {
    alpha1 <- support[1]
    alpha2 <- support[2]
    C1alf1 <- rho(Rinverse(alpha1))/alpha1
    Cbar1alf2 <- -rho(Rinverse(alpha2))/(1-alpha2)
    scorefn <- function(p, param) {
      stopifnot(p>=0, p<=1)
      if (param==1) {
        nu_low <- -C1alf1
        nu_high <- -Cbar1alf2
        mid_fn <- function(alf) lambdarho(Rinverse(alf))
      } else {
        nu_low <- -Rinverse(alpha1)*C1alf1
        nu_high <- -Rinverse(alpha2)*Cbar1alf2
        mid_fn <- function(p) Rinverse(p)*lambdarho(Rinverse(p))-1
      }
      Pregion <- cut(p, breaks=c(0, alpha1, alpha2, 1.01),
                     labels=c('low','mid','high'), include.lowest = TRUE)
      switch(as.character(Pregion),
             low = nu_low,
             mid = mid_fn(p),
             high = nu_high )
    }
    W <- list(purrr::map_dbl(PIT,~scorefn(.x,1)),
              purrr::map_dbl(PIT,~scorefn(.x,2)))
    return(W)
  }
  return(nu)
}


#' Probitnormal score kernel function
#'
#' Returns a function nu(P), where \eqn{\nu} itself returns a list with the
#' centered PNS kernels on the desired \code{support}.
#'
#' @inheritParams .tlsfkernel_function
#' @return \eqn{\nu(P)} kernel function
#' @export
nu_probitnormal <- function(support=c(0,1), param=NULL) {
  stopifnot(is.null(param))
  tlsf_scorefunction(support, qnorm, dnorm, identity)
}

# nu_probitnormal <- function(support=c(0,1), param=NULL) {
#   stopifnot(is.null(param))
#   nu <- function(PIT) {
#     alpha1 <- support[1]
#     alpha2 <- support[2]
#     scorefn <- function(p, param) {
#       stopifnot(p>=0, p<=1)
#       if (param==1) {
#         nu_low <- -dnorm(qnorm(alpha1))/alpha1
#         nu_high <- dnorm(qnorm(alpha2))/(1-alpha2)
#         mid_fn <- qnorm
#       } else {
#         nu_low <- -qnorm(alpha1)*dnorm(qnorm(alpha1))/alpha1
#         nu_high <- qnorm(alpha2)*dnorm(qnorm(alpha2))/(1-alpha2)
#         mid_fn <- function(p) qnorm(p)^2-1
#       }
#       Pregion <- cut(p, breaks=c(0, alpha1, alpha2, 1.01),
#                      labels=c('low','mid','high'), include.lowest = TRUE)
#       switch(as.character(Pregion),
#              low = nu_low,
#              mid = mid_fn(p),
#              high = nu_high )
#     }
#     W <- list(purrr::map_dbl(PIT,~scorefn(.x,1)),
#               purrr::map_dbl(PIT,~scorefn(.x,2)))
#     return(W)
#   }
#   return(nu)
# }


#' Probitnormal variance-covariance function
#'
#' @inheritParams .tlsfkernel_function
#' @return \eqn{\Sigma} variance matrix (2x2)
#' @export
vcv_probitnormal <- function(support=c(0,1), param=NULL) {
  stopifnot(is.null(param))
  Rinverse <- qnorm
  rho <- dnorm
  lambdarho <- identity
  # B integrals unique to each location-scale family
  B0 <- identity
  B1 <- function(alf) -rho(Rinverse(alf))
  B2 <- function(alf) if (alf==1) { 1 } else {
           alf+Rinverse(alf)*B1(alf) }
  vcv_tlsf(support, Rinverse, rho, lambdarho, B0, B1, B2)
}

#' Logit-Logistic score kernel function
#'
#' Returns a function nu(P), where \eqn{\nu} itself returns a list with the
#' centered LLS kernels on the desired \code{support}.
#'
#' @inheritParams .tlsfkernel_function
#' @return \eqn{\nu(P)} kernel function
#' @export
nu_logitlogistic <- function(support=c(0,1), param=NULL) {
  stopifnot(is.null(param))
  lambdarho <- function(x) 2*plogis(x)-1
  tlsf_scorefunction(support, qlogis, dlogis, lambdarho)
}

#' Logit-Logistic variance-covariance function
#'
#' @inheritParams .tlsfkernel_function
#' @return \eqn{\Sigma} variance matrix (2x2)
#' @export
vcv_logitlogistic <- function(support=c(0,1), param=NULL) {
  stopifnot(is.null(param))
  Rinverse <- qlogis
  rho <- dlogis
  lambdarho <- function(x) 2*plogis(x)-1
  # B integrals unique to each location-scale family
  B0 <- function(alf) alf^2*(3-2*alf)/3
  B1 <- function(alf) if (alf==1) { 0 } else {
    (alf*(1-alf)+log(1-alf))/3 + B0(alf)*Rinverse(alf) }
  # Asymptotic form from DLMF 25.12.3
  B2 <- function(alf) if (alf==1) {
    (2/3)*(pi^2/6-1)
    } else {
    (-2/3)*(alf+gsl::dilog(-alf/(1-alf))) + 2*B1(alf)*Rinverse(alf) -
     B0(alf)*Rinverse(alf)^2 }
  vcv_tlsf(support, Rinverse, rho, lambdarho, B0, B1, B2)
}

#' Logistic-beta score kernel function
#'
#' Returns a function nu(P; a,b), where \eqn{\nu} itself returns a list with the
#' centered LLS kernels on the desired \code{support}.
#'
#' @param support lower and upper bound of kernel window
#' @param param beta parameters (a,b)
#' @return \eqn{\nu(P)} kernel function
#' @export
nu_logisticbeta <- function(support=c(0,1), param=c(1,1)) {
  stopifnot(all(param>0))
  a<-param[1]
  b<-param[2]
#  Rinverse <- function(p) qlogis(qbeta(p,a,b))
  Rinverse <- function(p) {
    if (p<0.5) {
      qlogis(qbeta(p,a,b))
    } else {
      -qlogis(qbeta(1-p,b,a))
    }
  }
  rho <- function(x)
         plogis(x)^a*plogis(-x)^b/beta(a,b)
  lambdarho <- function(x) b*plogis(x)-a*plogis(-x)
  tlsf_scorefunction(support, Rinverse, rho, lambdarho)
}

#' Logistic-beta variance-covariance function
#'
#' @param support lower and upper bound of kernel window
#' @param param beta parameters (a,b)
#' @return \eqn{\Sigma} variance matrix (2x2)
#' @export
vcv_logisticbeta <- function(support=c(0,1), param=c(1,1)) {
  a<-max(param)
  b<-min(param)
  stopifnot(b>0)
  if (all(param==1)) return(vcv_logitlogistic(support))
  complementary <- (param[1]<param[2]) # flag
  R <- function(x) pbeta(plogis(x),a,b)
  Rinverse <- function(p) {
    if (p<0.5) {
      qlogis(qbeta(p,a,b))
    } else {
      -qlogis(qbeta(1-p,b,a))
    }
  }
  rho <- function(x)
    plogis(x)^a*plogis(-x)^b/beta(a,b)
  lambdarho <- function(x) b*plogis(x)-a*plogis(-x)
  # B integrals unique to each location-scale family
  B0alf1 <- a*b/(1+a+b)
  # kappa1 <- gsl::psi(a+1)-gsl::psi(b+1)
  # kappa2 <- gsl::psi_1(a+1)+gsl::psi_1(b+1)
  B0 <- function(alf)
    B0alf1*pbeta(plogis(Rinverse(alf)),a+1,b+1)
  if (a==3/2 & b==1/2) {
    B1B2 <- logisticBeta_3halves_1half
  } else if (a==1 & b==1/2) {
    B1B2 <- logisticBeta_one_1half
  } else if (a==2/3 & b==1/3) {
    B1B2 <- logisticBeta_2thirds_1third
  } else if (a==1 & b==1) {
    B1B2 <- logisticBeta_one_one
  } else if (a==1.02 & b==1.02) {
    B1B2 <- lBeta4
  } else if (a==0.98 & b==0.98) {
    B1B2 <- lBeta3
  } else if (a==3/2 & b==3/2) {
    B1B2 <- logisticBeta_3halves_3halves
  } else if (a==1/2 & b==1/2) {
    B1B2 <- logisticBeta_1half_1half
  } else if (a==1.02 & b==0.98) {
    B1B2 <- lBeta1
  } else {
    stop("Unrecognized LogisticBeta case")
  }
  B1 <- function(alf) with(B1B2, approx(alpha,B1,alf)$y)
  B2 <- function(alf) with(B1B2, approx(alpha,B2,alf)$y)
  if (complementary) {
    Rlist <- tlsf_complementary(R,Rinverse,rho,lambdarho, B0, B1, B2)
    vcv_tlsf(support, Rlist$Rinverse, Rlist$rho, Rlist$lambdarho,
             Rlist$B0, Rlist$B1, Rlist$B2)
  } else {
    vcv_tlsf(support, Rinverse, rho, lambdarho, B0, B1, B2)
  }
}

# Note that (skewness,kurtosis) close to Gumbel for (a,b)=(3/2,1/2).
# For (a,b)=(2/3,1/3), (skewness,kurtosis)=(1,6).

#' Cumulants of the Logistic-beta kernel on unit interval window
#'
#' @param support lower and upper bound of kernel window
#' @param param beta parameters (a,b)
#' @return list with elements \code{cumulants} (list length 4) and scalars \code{skewness} and \code{kurtosis}.
#' @export
cumulants_logisticbeta <- function(param=c(1,1)) {
  a<-param[1]
  b<-param[2]
  cumulants <- c(gsl::psi(a)-gsl::psi(b),
                 gsl::psi_1(a)+gsl::psi_1(b),
                 gsl::psi_n(2,a)-gsl::psi_n(2,b),
                 gsl::psi_n(3,a)+gsl::psi_n(3,b))

  skewness <- cumulants[3]/cumulants[2]^(3/2)
  kurtosis <- 3+cumulants[4]/cumulants[2]^2

  list(cumulants=cumulants, skewness=skewness, kurtosis=kurtosis)
}

#' Complementary distribution functions
#'
#' @param R cdf of distribution
#' @param Rinverse inverse cdf of distribution
#' @param rho density of distribution
#' @param lambdarho \eqn{\lambda_\rho} function
#' @param B0,B1,B2 integrals expressed as functions of alpha
#' @return named list of complementary functions
#' @export
tlsf_complementary <- function(R, Rinverse, rho, lambdarho,
                               B0=NULL, B1=NULL, B2=NULL) {
  if (!is.null(B0)) {
    B0c <- function(alf) B0(1)-B0(1-alf)
    B1c <- function(alf) -(B1(1)-B1(1-alf))
    B2c <- function(alf) B2(1)-B2(1-alf)
  } else {
    B0c <- NULL
    B1c <- NULL
    B2c <- NULL
  }

  list(
    R=function(x) 1-R(-x),
    Rinverse=function(p) -Rinverse(1-p),
    rho=function(x) rho(-x),
    lambdarho= function(x) -lambdarho(-x),
    B0=B0c, B1=B1c, B2=B2c)
}

#' loglog-Gumbel score kernel function
#'
#' Returns a function nu(P), where \eqn{\nu} itself returns a list with the
#' centered LLS kernels on the desired \code{support}.
#'
#' @param support lower and upper bound of kernel window
#' @param param (boolean) if TRUE, use complementary distribution
#' @return \eqn{\nu(P)} kernel function
#' @export
nu_gumbel <- function(support=c(0,1), param=FALSE) {
  R <- function(x) exp(-exp(-x))
  Rinverse <- function(p) -log(-log(p))
  rho <- function(x) exp(-x)*R(x)
  lambdarho <- function(x) 1-exp(-x)
  if (param) {
    Rlist <- tlsf_complementary(R,Rinverse,rho,lambdarho)
    tlsf_scorefunction(support, Rlist$Rinverse, Rlist$rho, Rlist$lambdarho)
  } else {
    tlsf_scorefunction(support, Rinverse, rho, lambdarho)
  }
}

#' loglog-Gumbel variance-covariance function
#'
#' @param support lower and upper bound of kernel window
#' @param param (boolean) if TRUE, use complementary distribution
#' @return \eqn{\Sigma} variance matrix (2x2)
#' @export
vcv_gumbel <- function(support=c(0,1), param=FALSE) {
  R <- function(x) exp(-exp(-x))
  Rinverse <- function(p) -log(-log(p))
  rho <- function(x) exp(-x)*R(x)
  lambdarho <- function(x) 1-exp(-x)
  # B integrals unique to each location-scale family
  B0 <- function(alf) if (alf==0) { 0 } else {alf*(1-log(alf)) }
  B1 <- function(alf) {
    if (alf==0) { 0 } else {
      if (alf==1) { -(1+digamma(1)) } else {
        gsl::expint_Ei(log(alf)) - alf + B0(alf)*Rinverse(alf) }
    }
  }
  # Asymptotic form from DLMF 25.12.3
  B2 <- function(alf) with(gumbelB2, approx(alpha,B2,alf)$y)
  if (param) {
    Rlist <- tlsf_complementary(R,Rinverse,rho,lambdarho, B0, B1, B2)
    vcv_tlsf(support, Rlist$Rinverse, Rlist$rho, Rlist$lambdarho,
             Rlist$B0, Rlist$B1, Rlist$B2)
  } else {
    vcv_tlsf(support, Rinverse, rho, lambdarho, B0, B1, B2)
  }
}

#
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

# #' Probitnormal score kernel1 moments
# #'
# #' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
# #' PNS kernel1 on the desired \code{support}.
# #'
# #' @inheritParams .monokernel_moments
# #' @export
# mu_probitnormal1 <- function(support=c(0,1), param=NULL) {
#   qq <- qnorm(support)
#   dd <- dnorm(qq)
#   mu2 <- (dd[1]^2)/support[1] + (dd[2]^2)/(1-support[2]) +
#     dd[1]*qq[1] - dd[2]*qq[2] + diff(support)
#   mu <- c(mu1=0, mu2=mu2)
#   return(mu)
# }
#
# #' Probitnormal score kernel2 moments
# #'
# #' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
# #' PNS kernel2 on the desired \code{support}.
# #'
# #' @inheritParams .monokernel_moments
# #' @export
# mu_probitnormal2 <- function(support=c(0,1), param=NULL) {
#   qq <- qnorm(support)
#   dd <- dnorm(qq)
#   mu2 <- dd[1]*qq[1] - dd[2]*qq[2] + dd[1]*qq[1]^3 - dd[2]*qq[2]^3 +
#     ((dd[1]*qq[1])^2)/support[1] + ((dd[2]*qq[2])^2)/(1-support[2]) +
#     2*diff(support)
#   mu <- c(mu1=0, mu2=mu2)
#   return(mu)
# }

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

# rho_probitnormal <- function(support=c(0,1),param=NULL) {
#   qq <- qnorm(support)
#   dd <- dnorm(qq)
#   mux <- dd[1]*(1+qq[1]^2) - dd[2]*(1+qq[2]^2) +
#     (qq[1]*dd[1]^2)/support[1] + (qq[2]*dd[2]^2)/(1-support[2])
#   rho <-cross_moment_to_correlation(mux,
#                                     mu_probitnormal1(support),
#                                     mu_probitnormal2(support))
#   return(rho)
# }
