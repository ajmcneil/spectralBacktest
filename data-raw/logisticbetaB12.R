# Interpolation table for LogisticBeta B1 and B2 functions

# Note: B1 is minimized at alpha=R(0,a,b).
B1B2 <- function(a,b, n=10000) {
  Q <- function(p) qlogis(qbeta(p,a,b))
  lambdaprimeQ <- function(p)
     (a+b)*qbeta(p, a,b)*qbeta(1-p,b,a)
  # B end-points for B1 and B2 integrals
  B0alf1 <- a*b/(1+a+b)
  kappa1 <- gsl::psi(a+1)-gsl::psi(b+1)
  kappa2 <- gsl::psi_1(a+1)+gsl::psi_1(b+1)
  B1alf1 <- B0alf1*kappa1
  B2alf1 <- B0alf1*(kappa2+kappa1^2)
  PIT <- seq(1/n,1-1/n,by=1/n)
  f1 <- Q(PIT)*lambdaprimeQ(PIT)
  f2 <- Q(PIT)*f1
  PIT <- c(0,PIT,1)
  y1 <- cumsum(c(0, f1/n, 0))
#  y1 <- y*B1alf1/y[length(y)]
  y <- cumsum(c(0, f2/n, 0))
  y2 <- y*B2alf1/y[length(y)]
  data.frame(alpha=PIT, B1=y1, B2=y2)
}

# Now use the function on desired (a,b) pairs
logisticBeta_3halves_1half <- B1B2(3/2,1/2)
usethis::use_data(logisticBeta_3halves_1half, overwrite=TRUE)
logisticBeta_one_1half <- B1B2(1,1/2)
usethis::use_data(logisticBeta_one_1half, overwrite = TRUE)
logisticBeta_2thirds_1third <- B1B2(2/3,1/3)
usethis::use_data(logisticBeta_2thirds_1third, overwrite = TRUE)

# Testing. Delete later.
logisticBeta_1half_3halves <- B1B2(1/2,3/2)
usethis::use_data(logisticBeta_1half_3halves, overwrite=TRUE)
logisticBeta_1third_2thirds <- B1B2(1/3,2/3)
usethis::use_data(logisticBeta_1third_2thirds, overwrite = TRUE)
lBeta1 <- B1B2(1.02,0.98)
usethis::use_data(lBeta1, overwrite=TRUE)
lBeta2 <- B1B2(0.98,1.02)
usethis::use_data(lBeta2, overwrite=TRUE)
logisticBeta_one_one <- B1B2(1,1)
usethis::use_data(logisticBeta_one_one, overwrite=TRUE)
logisticBeta_1half_1half <- B1B2(1/2,1/2)
usethis::use_data(logisticBeta_1half_1half, overwrite=TRUE)
logisticBeta_3halves_3halves <- B1B2(3/2,3/2)
usethis::use_data(logisticBeta_3halves_3halves, overwrite=TRUE)
lBeta3 <- B1B2(0.98, 0.98)
usethis::use_data(lBeta3, overwrite=TRUE)
lBeta4 <- B1B2(1.02, 1.02)
usethis::use_data(lBeta4, overwrite=TRUE)


