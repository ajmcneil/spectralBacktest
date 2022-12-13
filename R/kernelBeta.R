# Extended beta kernel

# Helper function for positive integer test
is.natural <- function(z)
  (z==round(z)) & (z>0)


# Helper function for brute-force series expansion of any hypergeometric
hypergeomSeries <- function(a,b,z,tol=1e-8) {
  # f <- pFq <- 1
  # # force j to be large enough to eliminate sign flipping
  # jstopmin <- min(a,b)+3
  # j <- 0
  # while ((abs(f/pFq) > tol) | (j<jstopmin)) {
  #   j <- j+1
  #   f <- f*z*prod(a+j-1)/(j*prod(b+j-1))
  #   pFq <- pFq+f
  # }
  # pFq
  Re(hypergeo::genhypergeo(a,b,z,tol=tol,maxiter=25000,series=TRUE))
}

# Helper function for the beta kernel.  Provides series solution
# based on
# M = (beta(a1+a2,1+b1+b2)/a2)*3F2([1, a1+a2+k, a2+b2], [1+a2, 1+a1+a2+b1+b2]; z)
# slow if a2 near 1
Mseries3F2 <- function(a1,b1,a2,b2,z=1,k=0,tol=1e-8) {
  # f <- F <- 1
  # j <- 0
  # while (f/F > tol) {
  #   j <- j+1
  #   f <- f*z*(a1+a2+k+j-1)*(a2+b2+j-1)/((a2+j)*(a1+a2+b1+b2+j))
  #   F <- F+f
  # }
  hypergeomSeries(c(1,a1+a2,a2+b2),
                  c(1+a2,1+a1+a2+b1+b2),
                  z, tol=tol)*beta(a1+a2,1+b1+b2)/a2
  # Re(hypergeo::genhypergeo(c(1,a1+a2,a2+b2),
  #                       c(1+a2,1+a1+a2+b1+b2),
  #                       z, tol=1e-12))*beta(a1+a2,1+b1+b2)/a2
}

# Note on gsl:
#  gsl_sf_beta() is a dropin for R beta().
#  gsl::beta_inc(a,b,x) is the same as pbeta() but handles b<0

# Helper function for incomplete beta when b==0 and 2*a is an integer.
beta_inc_B_0 <- function(x,a) {
  stopifnot("parameter a is integer or half-integer"=is.natural(2*a))
  halfint <- a %% 1
  h <- -log(1-x)
  if (halfint>0)
    h <- h + 2*log(1+sqrt(x))
  if (a>1.4) {
    w <- (1:ceiling(a-1))-halfint
    h <- h - drop(crossprod(1/w, x^w))
  }
  return(h)
}

#' Incomplete Beta B(x, a, b)
#'
#' Returns the incomplete beta B(x, a, b) integral.
#'
#' For parameters (a,b) both positive, this is equal to
#' pbeta(x,a,b)*beta(a,b). When b nonpositive, we use
#' the 2F1 function via DLMF 17.07.  We could alternatively
#' use gsl::beta_inc(a,b,x)*gsl::gsl_sf_beta(a,b), but the
#' gsl::hyperg_2F1() function is faster and handles b=0.
#' Unfortunately, for x near 1, the hyperg_2F1 solution
#' appears to be unstable at b==0.
#' @param x (vector) real number in unit interval
#' @param a (scalar) first parameter must be positive real.
#' @param b (scalar) second parameter must be real
#' @return B (vector) of incomplete beta evaluated at x
#' @export
beta_inc_B <- function(x,a,b) {
  stopifnot("argument must be in unit interval"= all(x>=0 & x<=1))
  stopifnot("parameter a must be positive"= (a>0))
  if (b>0) {
    pbeta(x,a,b)*beta(a,b)
  } else {
    if ((b==0) && is.natural(2*a)) {
      beta_inc_B_0(x,a)
    } else {
      # DLMF 17.07
      (x^a/a)*Re(hypergeo::hypergeo(a,1-b,a+1,x))
    }
  }
}

# Helper function: integral of (1-u)*g_1(u)*G_2(u)
mustar_beta <- function(p1, q1, p2=p1, q2=q1, z=1, force_Mseries=FALSE) {
  if (z==0) return(0)
  if (force_Mseries) return(Mseries3F2(p1,q1,p2,q2,z))
  if (is.natural(q2)) { # A.5: okay for z==1
    M <- 0
    for (k in 0:round(q2-1))
      M <- M + exp(gsl::lnpoch(p2,k)-gsl::lnfact(k))*beta_inc_B(z, p1+p2, 1+q1+k)
    M <- M*exp(gsl::lnfact(q2-1)-gsl::lnpoch(p2,q2))
  } else if (q1==0) {
    if ((z==1) & (q2==0)) {
      M <- (digamma(p1+p2)-digamma(p2))/p1 # B.1
    } else {  # A.2
      M <- (z^p1*beta_inc_B(z,p2,q2)-beta_inc_B(z,p1+p2,q2))/p1
    }
  } else if (is.natural(q1)) { # A.3 which relies on A.2, B.1
    M <- 0
    for (k in 0:round(q1))
      M <- M + (-1)^k*choose(q1,k)*mustar_beta(p1+k,0,p2,q2,z)
  } else if (is.natural(p2)) { # A.4. z==1 is okay
    if (p2==1) {
      if (q2==0) {
        # Initialize to the value when z==1
        M <- beta(p1,1+q1)*(digamma(1+p1+q1)-digamma(1+q1))
        if (z<1) {
          M <- M + beta_inc_B(1-z,1+q1,p1)*log(1-z) -
            (1/(1+q1)^2)*hypergeomSeries(c(1-p1,1+q1,1+q1),
                                         c(2+q1,2+q1),
                                         1-z, tol=1e-12)*(1-z)^(1+q1)
        }
      } else {  #q2 neq 0
        M <- (beta_inc_B(z,p1,1+q1)-beta_inc_B(z,p1,1+q1+q2))/q2
      }
    } else { # recursion for p2>1
      # M <- ((p2-1)/(p2+q2-1))*(mustar_beta(p1,q1,p2-1,q2,z)-
      #                       beta_inc_B(z,p1+p2-1,1+q1+q2))
      M <- ((p2-1)*mustar_beta(p1,q1,p2-1,q2,z)-
              beta_inc_B(z,p1+p2-1,1+q1+q2))/(p2+q2-1)
    }
  } else if (z<1) {
    m <- q1+q2
    if (m==round(m)) { # A.1
      M <- 0
      f <- 1
      for (k in 0:m) {
        M <- M + f*hypergeomSeries(c(1,p1+p2+k,p2+q2),
                                   c(1+p2, 1+p1+p2+m),z)
        f <- f*(1-z)*(p1+p2+k)/(k+1)
      }
      M <- M*z^(p1+p2)*beta(p1+p2,1+m)/p2
    } else {
      M <- NaN  # got nothing for this case!
    }
  } else {  # z==1
    M <- Mseries3F2(p1,q1,p2,q2)
  }
  return(M)
}

#' Beta kernel function
#'
#' Returns a function nu(P), where nu is the beta kernel on
#' the desired \code{support} with beta parameters \eqn{(p,q)=}\code{param$beta}
#' and scaling parameters \code{param$scaling}.  If \code{param} is a vector,
#' we set \code{scaling} to the \code{support}.
#'
#' @inheritParams .monokernel_function
#' @export
nu_beta <- function(support=c(0,1), param, standardize=TRUE) {
  nu <- function(PIT) {
    if (is.list(param)) {
      p <- param$beta[1]
      q <- param$beta[2]
      scaling <- param$scaling
    } else {
      p <- param[1]
      q <- param[2]
      scaling <- support
    }
    stopifnot(0<=scaling[1] & scaling[1]<=support[1],
              1>=scaling[2] & scaling[2]>=support[2],
              q>0 | scaling[2]==1 | support[2]<scaling[2])
    # Rescale support
    scaledsupport <- (support-scaling[1])/diff(scaling)
    PITscaled <- truncatePIT(PIT, scaling, rescale=TRUE)
    W <- beta_inc_B(pmin(PITscaled, scaledsupport[2]), p,q)-
      beta_inc_B(pmin(PITscaled, scaledsupport[1]), p,q)
    if (standardize) {
      mu <- mu_beta(support, c(p,q), scaling) %>% unname()
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}
#' Beta moments
#'
#' Returns a vector of the first two uncentered moments \eqn{(\mu_1,\mu_2)} for the
#' Beta kernel scaled to \code{scaling} on the desired \code{support}
#' with parameter \eqn{(p,q)=}\code{param}.  Observe the divergence from
#' nu_beta() in calling syntax.
#'
#' @inheritParams .monokernel_moments
#' @param scaling lower and upper bounds on which to rescale the beta
#' @export
mu_beta <- function(support=c(0,1), param=c(1,1), scaling=support,
                    force_Mseries=FALSE) {
  p <- param[1]
  q <- param[2]
  scaledsupport <- (support-scaling[1])/diff(scaling)
  Gtm <- function(dq)
    diff(beta_inc_B(scaledsupport,p,q+dq))
  mu1 <- Gtm(1)
  mu2 <- 2*(mustar_beta(p,q,p,q,scaledsupport[2],force_Mseries)-
              mustar_beta(p,q,p,q,scaledsupport[1],force_Mseries)-
              beta_inc_B(scaledsupport[1],p,q)*mu1)
  if (scaling[2]==1) {
    totalmass <- 0  # avoid Inf*0 undefined
  } else {
    totalmass <- Gtm(0)
  }
  mu <- (1-scaling[2])*totalmass*c(1,totalmass) +
    diff(scaling)*c(mu1=mu1, mu2=mu2)
  return(mu)
}

#' Beta/Beta correlation
#'
#' Returns the correlation for the Beta/Beta pair on common \code{support}.
#'
#' @inheritParams .bikernel_correlation
#' @export
rho_beta_beta <- function(support=c(0,1), param, scaling=support) {
  p1 <- param[[1]][1]
  q1 <- param[[1]][2]
  p2 <- param[[2]][1]
  q2 <- param[[2]][2]
  scaledsupport <- (support-scaling[1])/diff(scaling)
  muA <- diff(beta_inc_B(scaledsupport,p1,q1+1))
  muB <- diff(beta_inc_B(scaledsupport,p2,q2+1))
  mux0 <- mustar_beta(p1,q1,p2,q2, scaledsupport[2]) -
    mustar_beta(p1,q1,p2,q2, scaledsupport[1]) +
    mustar_beta(p2,q2,p1,q1, scaledsupport[2]) -
    mustar_beta(p2,q2,p1,q1, scaledsupport[1]) -
    (muA*beta_inc_B(scaledsupport[1],p2,q2) +
       muB*beta_inc_B(scaledsupport[1],p1,q1))
  if (scaling[2]==1) {
    crossmass <- 0  # avoid Inf*0 undefined
  } else {
    crossmass <- diff(beta_inc_B(scaledsupport,p1,q1))*
      diff(beta_inc_B(scaledsupport,p2,q2))
  }
  cross_moment_to_correlation((1-scaling[2])*crossmass+diff(scaling)*mux0,
                              mu_beta(support, param[[1]], scaling),
                              mu_beta(support, param[[2]], scaling))
}

