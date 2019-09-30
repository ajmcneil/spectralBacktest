#' Generic singleton kernel function for inheritance of parameter documentation
#'
#' @param support lower and upper bound of kernel window
#' @param param optional parameters to kernel
#' @param standardize (boolean) if TRUE, then standardize \eqn{W} to mean 0, variance 1 under the null
#' @return \eqn{\nu(P)} kernel function
.monokernel_function <- function(support, param, standardize=TRUE)
  NULL

#' Generic singleton kernel moment function for inheritance of parameter documentation
#'
#' @param support lower and upper bound of kernel window
#' @param param optional parameters to kernel
#' @return vector of first two uncentered moments
.monokernel_moments <- function(support, param)
  NULL

# #' Generic bikernel cross-moment function for inheritance of parameter documentation
# #'
# #' @param support lower and upper bound of kernel window
# #' @param param list of optional parameters to kernels
# #' @return cross-moment
# .bikernel_crossmoment <- function(support, param)
#  NULL

#' Generic bikernel correlation function for inheritance of parameter documentation
#'
#' @param support lower and upper bound of kernel window
#' @param param list of optional parameters to kernels
#' @return correlation
.bikernel_correlation <- function(support, param)
  NULL

#' Truncate PIT to a kernel window and (optionally) rescale to unit interval.
#'
#' Calculate vector \eqn{P_t^* = \alpha_1 \vee (P_t \wedge \alpha_2)}.
#' If rescale==TRUE, then rescale to unit interval.
#'
#' @param PIT Vector of probability integral transform values
#' @param support lower and upper bound of kernel window
#' @return PITstar is vector of truncated PIT.
#' @export
truncatePIT <- function(PIT,support=c(0,1),rescale=FALSE) {
  Pstar <- pmax(pmin(PIT,support[2]),support[1])
  if (rescale)
    Pstar <- (Pstar-support[1])/diff(support)
  return(Pstar)
}

#' Epanechnikov kernel function
#'
#' Returns a function nu(P), where nu is the Epanechnikov kernel on
#' the desired \code{support}.
#'
#' @inheritParams .monokernel_function
#' @export
nu_epanechnikov <- function(support=c(0,1), param=NULL, standardize=TRUE) {
  nu <- function(PIT) {
    P <- truncatePIT(PIT, support, rescale = TRUE)
    W <- P^2*(3-2*P)
    if (standardize) {
      mu <- as.numeric(mu_epanechnikov(support))
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}

#' Epanechnikov moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' Epanechnikov kernel on the desired \code{support}.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_epanechnikov <- function(support=c(0,1),param=NULL) {
  1-support[2]+diff(support)*c(mu1=1/2, mu2=13/35)
}

#' Uniform kernel function
#'
#' Returns a function nu(P), where nu is the uniform kernel on
#' the desired \code{support}.
#'
#' @inheritParams .monokernel_function
#' @export
nu_uniform <- function(support=c(0,1), param=NULL, standardize=TRUE) {
  nu <- function(PIT) {
    W <- truncatePIT(PIT, support, rescale=TRUE)
    if (standardize) {
      mu <- as.numeric(mu_uniform(support))
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}

#' Uniform moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' uniform kernel on the desired \code{support}.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_uniform <- function(support=c(0,1),param=NULL) {
  1-support[2]+diff(support)*c(mu1=1/2, mu2=1/3)
}


#' Linear kernel function
#'
#' Returns a function nu(P), where nu is the linear kernel on
#' the desired \code{support}. The slope of the kernel density is given
#' by \code{param}, which takes values in {-1,1}.
#'
#' @inheritParams .monokernel_function
#' @export
nu_linear <- function(support=c(0,1), param=1, standardize=TRUE) {
  nu <- function(PIT) {
    P <- truncatePIT(PIT, support, rescale=TRUE)
    W <- P*(1-sign(param)*(1-P))
    if (standardize) {
      mu <- as.numeric(mu_linear(support,param))
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}

#' Linear moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' Linear kernel on the desired \code{support}.The slope of the kernel density is given
#' by \code{param}, which takes values in {-1,1}.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_linear <- function(support=c(0,1),param=1) {
  if (param>0) {
    mu <- c(mu1=1/3, mu2=1/5)
  } else{
    mu <- c(mu1=2/3, mu2=8/15)
  }
  1-support[2]+diff(support)*mu
}

#' Arcsin kernel function
#'
#' Returns a function nu(P), where nu is the arcsin kernel on
#' the desired \code{support}. Equivalent to nu_beta(param=c(1/2,1/2)).
#'
#' @inheritParams .monokernel_function
#' @export
nu_arcsin <- function(support=c(0,1), param=NULL, standardize=TRUE) {
  nu <- function(PIT) {
    W <- pbeta(truncatePIT(PIT, support, rescale=TRUE),1/2,1/2)
    if (standardize) {
      mu <- as.numeric(mu_arcsin(support,param))
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}

#' Arcsin moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' arcsin kernel on the desired \code{support}.Equivalent to mu_beta(param=c(1/2,1/2)).
#'
#' @inheritParams .monokernel_moments
#' @export
mu_arcsin <- function(support=c(0,1),param=NULL) {
  1-support[2]+diff(support)*c(1/2,1/2-2/pi^2)
}

#' Exponential kernel function
#'
#' Returns a function nu(P), where nu is the exponential kernel on
#' the desired \code{support} with parameter \eqn{\zeta=}\code{param}.
#'
#' @inheritParams .monokernel_function
#' @export
nu_exponential <- function(support=c(0,1), param=0, standardize=TRUE) {
  nu <- function(PIT) {
    if (param==0)
      W <- nu_uniform(support)
    else {
      P <- truncatePIT(PIT, support, rescale=TRUE)
      W <- (exp(param*P)-1)*diff(support)/param
    }
    if (standardize) {
      mu <- as.numeric(mu_exponential(support, param))
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}

#' Exponential moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' Exponential kernel on the desired \code{support} with parameter \eqn{\zeta=}\code{param}.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_exponential <- function(support=c(0,1), param=0) {
  if (param==0) {
    mu <- mu_uniform(support)
  } else {
    a <- support[1]
    b <- support[2]
    k <- param/(b-a)
    mu1 <- -(b*k - k - 1)*exp(-a*k + b*k)/k^2 + ((a - 1)*k - 1)/k^2
    mu2 <- -1/2*((2*b*k - 2*k - 1)*exp(2*b*k) - 4*(b*k*exp(a*k) -
                         (k + 1)*exp(a*k))*exp(b*k))*exp(-2*a*k)/k^3 - 1/2*(2*(a - 1)*k - 3)/k^3
    mu <- c(mu1=mu1, mu2=mu2)
  }
  return(mu)
}

# Helper function for the beta kernel.  Provides
# 3F2([1, a1+a2, a2+b2], [1+a2, 1+a1+a2+b1+b2]; 1)
series3F2 <- function(a1,b1,a2,b2,tol=1e-14) {
  f <- F <- 1
  k <- 0
  while (f/F > tol) {
    k <- k+1
    f <- f*(a1+a2+k-1)*(a2+b2+k-1)/((a2+k)*(a1+a2+b1+b2+k))
    F <- F+f
  }
  return(F)
}

#' Beta kernel function
#'
#' Returns a function nu(P), where nu is the beta kernel on
#' the desired \code{support} with parameter \eqn{(p,q)=}\code{param}.
#'
#' @inheritParams .monokernel_function
#' @export
nu_beta <- function(support=c(0,1), param=c(1,1), standardize=TRUE) {
  nu <- function(PIT) {
    W <- pbeta(truncatePIT(PIT, support, rescale=TRUE),param[1],param[2])
    if (standardize) {
      mu <- as.numeric(mu_beta(support, param))
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}

#' Beta moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' Beta kernel on the desired \code{support} with parameter \eqn{(p,q)=}\code{param}.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_beta <- function(support=c(0,1), param=c(1,1)) {
  p <- param[1]
  q <- param[2]
  mu1 <- q/(p+q)
  mu2 <- 2*(beta(2*p,1+2*q)/(p*beta(p,q)^2))*series3F2(p,q,p,q)
  mu <- 1-support[2]+diff(support)*c(mu1=mu1, mu2=mu2)
  return(mu)
}

#' Discrete kernel function
#'
#' Returns a function nu(P), where nu is the discrete kernel on
#' the desired \code{support} with weights in the \code{param} vector.
#'
#' @inheritParams .monokernel_function
#' @export
nu_discrete <- function(support=0.99, param=rep(1,length(support)), standardize=TRUE) {
  nu <- function(PIT) {
    violations <- mapply(function(alpha,wght) wght*exceedancePIT(PIT,alpha),
                         support, param)
    W <- rowSums(violations)
    if (standardize) {
      mu <- as.numeric(mu_discrete(support, param))
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}

#' Discrete kernel moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' discrete kernel on the desired \code{support} with weights in the \code{param} vector.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_discrete <- function(support=0.99, param=rep(1,length(support))) {
  mu1 <- param %*% (1-support)
  mu2 <- (2*param*cumsum(param)- param^2) %*% (1-support)
  mu <- c(mu1=as.numeric(mu1), mu2=as.numeric(mu2))
  return(mu)
}

#' Binomial kernel function for use in Pearson test.
#'
#' Returns a function nu(P), where nu is an indicator for PIT>\code{param}.
#' Note that \code{support} is ignored, and \code{param} takes its usual place.
#' This is to allow the Pearson test to be shoehorned into
#' the multikernel testing framework.
#'
#' @inheritParams .monokernel_function
#' @export
nu_pearson <- function(support=NULL, param=0.99, standardize=TRUE) {
  nu <- function(PIT) {
    W <- exceedancePIT(PIT,param)
    if (standardize) {
      mu <- as.numeric(mu_pearson(support, param))
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}

#' Binomial kernel moments for the Pearson test.
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' binomial kernel for PIT>\code{param}.
#'
#' Note that \code{support} is ignored, and \code{param} takes its usual place.
#' This is to allow the Pearson test to be shoehorned into the multikernel
#' testing framework.on the desired \code{support} with weights in the \code{param} vector.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_pearson <- function(support=NULL, param=0.99) {
  mu2 <- mu1 <- as.numeric(1-param)
  mu <- c(mu1=mu1, mu2=mu2)
  return(mu)
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
  as.numeric((mux - muA[1]*muB[1])/sigma)
}


#' Pearson correlation
#'
#' Returns a bivariate correlation for VaR exceedance pairs for use in the Pearson test.
#'
#' Note that \code{support} is ignored, and \code{param} takes its usual place.
#' This is to allow the Pearson test to be shoehorned into the multikernel
#' testing framework.on the desired \code{support} with weights in the \code{param} vector.
#'
#' @inheritParams .bikernel_correlation
#' @export
rho_pearson <- function(support=NULL, param=list(0.99,0.99)) {
  levels <- unlist(param)
  cross_moment_to_correlation(1-max(levels),
                              mu_pearson(param=levels[1]),
                              mu_pearson(param=levels[2]))
}

#' Linear/Linear correlation (positive/negative)
#'
#' Returns the correlation for the Linear/Linear pair on common
#' \code{support}.
#'
#' @inheritParams .bikernel_correlation
#' @export
rho_linear_linear <- function(support=c(0,1), param=NULL)
  cross_moment_to_correlation(1-support[2]+diff(support)*0.3,
                              mu_linear(support, param=-1),
                              mu_linear(support, param=1) )

#' Arcsin/Epanechnikov correlation
#'
#' Returns the correlation for the Arcsin/Epanechnikov pair on common
#' \code{support}.
#'
#' @inheritParams .bikernel_correlation
#' @export
rho_arcsin_epanechnikov <- function(support=c(0,1), param=NULL)
  cross_moment_to_correlation(1-support[2]+diff(support)*83/256,
                              mu_arcsin(support),
                              mu_epanechnikov(support) )


#' Beta/Beta correlation
#'
#' Returns the correlation for the Beta/Beta pair on common \code{support}.
#'
#' @inheritParams .bikernel_correlation
#' @export
rho_beta_beta <- function(support=c(0,1), param=list(NULL,NULL)) {
  alf1 <- support[1]
  alf2 <- support[2]
  p1 <- param[[1]][1]
  q1 <- param[[1]][2]
  p2 <- param[[2]][1]
  q2 <- param[[2]][2]
  betratio <- beta(p1+p2,1+q1+q2)/(beta(p1,q1)*beta(p2,q2))  # same for both kernels
  mux <- betratio*(series3F2(p1,q1,p2,q2)/p2 + series3F2(p2,q2,p1,q1)/p1)
  cross_moment_to_correlation(1-alf2+(alf2-alf1)*mux,
                              mu_beta(support,param[[1]]),
                              mu_beta(support,param[[2]]))
}

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

#' Probitnormal score kernel function
#'
#' Returns a function nu(P), where \eqn{\nu} is the centered PNS j-kernel (\eqn{j=1,2}) on
#' the desired \code{support}. Set \code{param} to kernel \eqn{j}.
#'
#' @inheritParams .monokernel_function
#' @export
nu_probitnormal <- function(support=c(0,1), param=0, standardize=TRUE) {
  nu <- function(PIT) {
    stopifnot(param %in% c(1,2))
    alpha1 <- support[1]
    alpha2 <- support[2]
    if (param==1) {
      nu_low <- -dnorm(qnorm(alpha1))/alpha1
      nu_high <- dnorm(qnorm(alpha2))/(1-alpha2)
      mid_fn <- qnorm
    } else {
      nu_low <- -qnorm(alpha1)*dnorm(qnorm(alpha1))/alpha1
      nu_high <- qnorm(alpha2)*dnorm(qnorm(alpha2))/(1-alpha2)
      mid_fn <- function(p) qnorm(p)^2-1
    }
    scorefn <- function(p) {
      Pregion <- cut(p, breaks=c(0, alpha1, alpha2, 1),
                     labels=c('low','mid','high'), include.lowest = TRUE)
      switch(as.character(Pregion),
             low = nu_low,
             mid = mid_fn(p),
             high = nu_high )
    }
    W <- sapply(PIT,scorefn)
    if (standardize) {
      mu <- mu_probitnormal(support, param)
      W <- as.numeric((W-mu[1])/sqrt(mu[2]-mu[1]^2))
    }
    return(W)
  }
  return(nu)
}

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

#' Probitnormal score kernel moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' PNS kernel on the desired \code{support}.  Set \code{param} to control kernel \eqn{j}.
#'
#' @inheritParams .monokernel_moments
#' @export
mu_probitnormal <- function(support=c(0,1), param=0) {
  stopifnot(param %in% c(1,2))
  qq <- qnorm(support)
  dd <- dnorm(qq)
  if (param==1) {
    mu2 <- (dd[1]^2)/support[1] + (dd[2]^2)/(1-support[2]) +
      dd[1]*qq[1] - dd[2]*qq[2] + diff(support)
  } else {
    mu2 <- dd[1]*qq[1] - dd[2]*qq[2] + dd[1]*qq[1]^3 - dd[2]*qq[2]^3 +
            ((dd[1]*qq[1])^2)/support[1] + ((dd[2]*qq[2])^2)/(1-support[2]) +
            2*diff(support)
  }
  mu <- c(mu1=0, mu2=mu2)
  return(mu)
}

#' Probitnormal score correlation
#'
#' Returns the correlation for VaR exceedance pairs for use in the PNS test.
#'
#' @inheritParams .bikernel_correlation
#' @export
rho_probitnormal <- function(support=c(0,1),param=NULL) {
  qq <- qnorm(support)
  dd <- dnorm(qq)
  mux <- dd[1]*(1+qq[1]^2) - dd[2]*(1+qq[2]^2) +
    (qq[1]*dd[1]^2)/support[1] + (qq[2]*dd[2]^2)/(1-support[2])
  rho <-cross_moment_to_correlation(mux,
                                    mu_probitnormal(support, param=1),
                                    mu_probitnormal(support, param=2))
  return(rho)
}

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

