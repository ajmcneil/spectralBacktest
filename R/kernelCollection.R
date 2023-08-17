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
      mu <- mu_epanechnikov(support) |> unname()
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
      mu <- mu_uniform(support) |> unname()
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
      mu <- mu_linear(support,param) |> unname()
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
      mu <- mu_arcsin(support,param) |> unname()
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
      mu <- mu_exponential(support, param) |> unname()
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
      mu <- mu_discrete(support, param) |> unname()
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
      mu <- mu_pearson(support, param) |> unname()
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
  mu2 <- mu1 <- 1-unname(param)
  mu <- c(mu1=mu1, mu2=mu2)
  return(mu)
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



