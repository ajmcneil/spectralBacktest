# Tests to confirm that the new extended beta moments are consistent with the old.
library(gsl)

# log Pochhammer function (n)_k.
lpochhammer <- function(n, k)
  lgamma(n+k)-lgamma(n)

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

# Old version of mustar_beta
mustar_beta_regularized <- function(p1, q1, p2, q2) {
  if (q2==1) {
    f <- q1*exp(lpochhammer(p1,p2)-lpochhammer(p1+q1,p2+1))
  } else if (p2==1) {
    f <- (q1/(p1+q1))-exp(lpochhammer(q1,q2+1)-lpochhammer(p1+q1,q2+1))
  } else {
    f <- (beta(p1+p2,1+q1+q2)/(p2*beta(p1,q1)*beta(p2,q2)))*series3F2(p1,q1,p2,q2)
  }
  f
}

# Old version of beta kernel
nu_beta_regularized <- function(support=c(0,1), param=c(1,1), standardize=TRUE) {
  nu <- function(PIT) {
    W <- pbeta(truncatePIT(PIT, support, rescale=TRUE),param[1],param[2])
    if (standardize) {
      mu <- mu_beta_regularized(support, param) %>% unname()
      W <- (W-mu[1])/sqrt(mu[2]-mu[1]^2)
    }
    return(W)
  }
  return(nu)
}

# Old version of beta kernel moments
mu_beta_regularized <- function(support=c(0,1), param=c(1,1)) {
  p <- param[1]
  q <- param[2]
  mu1 <- q/(p+q)
  mu2 <- 2*mustar_beta_regularized(p,q,p,q)
  mu <- 1-support[2]+diff(support)*c(mu1=mu1, mu2=mu2)
  return(mu)
}

# Old version deregularized manually
mu_beta_old <- function(support=c(0,1), param=c(1,1)) {
  bb <- beta(param[1],param[2])
  mu_beta_regularized(support, param)*c(bb,bb*bb)
}

# The first set of tests compare the old regularized beta kernel to the new kernel,
# which is not regularized. We need to consider p=1 and q=1 as special cases, and also
# the special case of alpha2=1.
tol <- 1e-5
suppt <- c(0.3,0.9)
test_that("Case of alpha2<1", {
  expect_equal(mu_beta(support=suppt, param=c(1,1)),
               mu_beta_old(support=suppt, param=c(1,1)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(2,1)),
               mu_beta_old(support=suppt, param=c(2,1)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(1.5,0.5)),
               mu_beta_old(support=suppt, param=c(1.5,0.5)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(1.5,1)),
               mu_beta_old(support=suppt, param=c(1.5,1)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(1,0.5)),
               mu_beta_old(support=suppt, param=c(1,0.5)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(2.4,1.3)),
               mu_beta_old(support=suppt, param=c(2.4,1.3)), tolerance=tol, scale=1)
})
suppt[2] <- 1
test_that("Case of alpha2==1", {
  expect_equal(mu_beta(support=suppt, param=c(1,1)),
               mu_beta_old(support=suppt, param=c(1,1)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(2,1)),
               mu_beta_old(support=suppt, param=c(2,1)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(1.5,0.5)),
               mu_beta_old(support=suppt, param=c(1.5,0.5)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(1.5,1)),
               mu_beta_old(support=suppt, param=c(1.5,1)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(1,0.5)),
               mu_beta_old(support=suppt, param=c(1,0.5)), tolerance=tol, scale=1)
  expect_equal(mu_beta(support=suppt, param=c(2.4,1.3)),
               mu_beta_old(support=suppt, param=c(2.4,1.3)), tolerance=tol, scale=1)
})

# Need a set of tests for rho_beta_beta against other special-case bikernels
# Key question is whether scaling alters the correlation.  Is it perhaps
# enough to know the canonical correlation for scaling=c(0,1)?
