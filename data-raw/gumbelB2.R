# Interpolation table for Gumbel B2 function

R <- function(x) exp(-exp(-x))
Q <- function(p) -log(-log(p))
eulergamma <- -digamma(1)
#alflower <- R(0)
#B2lower <- 0.634454
alflower <- B2lower <- 0
n <- 10000
PIT <- alflower + (1-alflower)*seq(1/n,1-1/n,by=1/n)
f <- -log(PIT)*Q(PIT)^2
PIT <- c(alflower,PIT,1)
y <- cumsum(c(B2lower, f*(1-alflower)/n, 0))
y1 <- pi^2/6+eulergamma^2-2*eulergamma
y <- y*y1/y[length(y)]
gumbelB2 <- data.frame(alpha=PIT, B2=y)
usethis::use_data(gumbelB2, overwrite=TRUE)

# Interpolation dataframe for B1 and B2 functions.
#
# @param R cdf of TLSF
# @param Rinverse quantile function of TLSF \eqn{R}
# @param val1 value of integral at alpha=1 (if known)
# @param k exponent in \eqn{B_k} integral
# @param n number of rows in dataframe
# @return \eqn{I(\theta)} information matrix
# function tlsf_interp_table(R,Rinverse,val1=NULL, k=2,n=10000) {
#
# }
