# Functions for martingale difference tests
#
# WLOG, CVT h(P) scaled to have mean 0, variance 1 when P uniform


#' One-sided upper-tail CVT
#'
#' h(P) scaled to have mean 0, variance 1 when P uniform
#'
#' @param param level for exceedance indicator
#' @return \eqn{h(P)} MD-lag function
#' @export
CVT_uppertail <- function(param=0.99) {
  h <- function(P) {
    (exceedancePIT(P,param) - (1-param))/sqrt(param*(1-param))
  }
  return(h)
}

#' Two-sided tail MD-lag function
#'
#' Indicator for abs(2 P-1) > abs(2*level-1)
#' h(P) scaled to have mean 0, variance 1 when P uniform
#'
#' @param param level for exceedance indicator
#' @return \eqn{h(P)} MD-lag function
#' @export
CVT_twotail <- function(param=0.99) {
  h <- function(P) {
    level <- abs(2*param-1)
    return((exceedancePIT(abs(2*P-1),level) - (1-level))/sqrt(level*(1-level)))
  }
  return(h)
}

#' Power MD-lag function
#'
#' h(P) scaled to have mean 0, variance 1 when P uniform
#'
#' @param param power to which \eqn{\abs(2P-1)} is raised
#' @return \eqn{h(P)} MD-lag function
#' @export
CVT_PtildePower <- function(param=4) {
  h <- function(P) {
    ptilde <- abs(2*P-1)
    mu <- 1/(c(1,2)*param+1)
    return((ptilde^param-mu[1])/sqrt(mu[2]-mu[1]^2))
  }
  return(h)
}

#' Correlation of power MD-lag functions
#'
#' @param param list of powers to which \eqn{|2P-1|} is raised
#' @return \eqn{\rho} Correlation of \eqn{h_1(P)} and \eqn{h_2(P)} under null.
#' @export
rho_PtildePower <- function(param) {
  mu <- 1/(unlist(param)+1)
  s2 <- 1/(2*unlist(param)+1)-mu^2
  mux <- 1/(sum(unlist(param))+1)
  rho <- (mux-prod(mu))/sqrt(prod(s2))
  return(rho)
}

#' Matrix of lagged values of a vector
#'
#' @param x vector
#' @param max.lag number of lags to include
#' @return matrix of size (length(x), max.lag+1) of 0-order to max.lag-order lags of x
lagmatrix <- function(x,max.lag)
   embed(c(rep(NA,max.lag),x),max.lag+1)[,-1]

#' Empirical H function from h_matrix or h_list input
#'
#' @param h_matrix matrix of lagged \eqn{h} values (including constant column) or list of same
#' @return estimated H matrix
empiricalH <- function(h_matrix) {
  if (is.list(h_matrix))
      h_matrix <- do.call(cbind,h_matrix)
  h <- h_matrix[!is.na(rowSums(h_matrix)),]
  (t(h) %*% h)/nrow(h)
}

# Helper function to build a block of the partial information H matrix
buildblockH <- function(v12,k1,m1,k2=k1,m2=m1)
  diag(c(0, rep(v12, min(k1,k2))),k1+1, k2+1) +
    (c(1, rep(m1,k1)) %o% c(1, rep(m2,k2)))


#' Partial information H function from h_matrix or h_list input
#'
#' @param h_matrix matrix of lagged \eqn{h} values (including constant column) or list of same
#' @return estimated H matrix
#' @export
partialinfoH <- function(h_matrix) {
  if (is.list(h_matrix)) {
      stopifnot(length(h_matrix)==2)  # only bikernel is implemented
      kk <- purrr::map_dbl(h_matrix, ~NCOL(.x)-1)
      mink1 <- min(kk)+1
      hcolmean <- purrr::map(h_matrix,~colMeans(.x, na.rm=TRUE))
      hcolvar <- purrr::map2(h_matrix,hcolmean,
                            ~(colMeans(.x^2, na.rm=TRUE) - .y^2))
      hmean <- purrr::map_dbl(hcolmean,~mean(.x[-1]))
      hvar <- purrr::map_dbl(hcolvar,~mean(.x[-1]))
      hmink1 <- purrr::map(h_matrix, ~as.matrix(.x[,1:mink1]))
      hcrosscol <- colMeans(hmink1[[1]]*hmink1[[2]], na.rm=TRUE) -
                    hcolmean[[1]][1:mink1]*hcolmean[[2]][1:mink1]
      hcovar <- mean(hcrosscol[-1])
      A11 <- buildblockH(hvar[1], kk[1], hmean[1])
      A22 <- buildblockH(hvar[2], kk[2], hmean[2])
      A12 <- buildblockH(hcovar, kk[1], hmean[1], kk[2], hmean[2])
      rbind(cbind(A11, A12), cbind(t(A12), A22))
  } else {
    hcolmean <- colMeans(h_matrix, na.rm=TRUE)
    hcolvar <- colMeans(h_matrix^2, na.rm=TRUE) - hcolmean^2
    buildblockH(mean(hcolvar[-1]), NCOL(h_matrix)-1,
           mean(hcolmean[-1]))
  }
}

#' Theoretical H function from CVT input for bispectral and multispectral MD tests.
#'
#' @param CVT list describing lag structure
#' @return theoretical H matrix
theoreticalH <- function(CVT) {
   if (length(CVT$lags)==1) {
     H <- diag(1+CVT$lags)
   } else {
     if (length(CVT$lags)==2) {
       h_rho <- CVT$correlation(CVT$h_param)
       minlag <- min(CVT$lags)
       B <- diag(c(1, rep(h_rho, minlag)),
                 nrow=(1+CVT$lags[1]), ncol=(1+CVT$lags[2]))
       H <- rbind(cbind(diag(1+CVT$lags[1]), B),
                  cbind(t(B),diag(1+CVT$lags[2])))
     }
   }
}

#' Monospectral h matrix
#'
#' @param CVT list describing lag structure
#' @param PIT vector of PIT values (length \eqn{n})
#' @return \eqn{h} matrix (size n x (1+monokernel_CVT$lags))
#' @export
monospectral_h_matrix <- function(CVT, PIT)
   cbind(rep(1,length(PIT)),
         lagmatrix(CVT$h_closure(CVT$h_param)(PIT), CVT$lags))


#' Multispectral h matrix
#'
#' @param CVT list describing bikernel or multikernel lag structure
#' @param PIT vector of PIT values (length \eqn{n})
#' @return list of \eqn{h} matrices (size n x (1+CVT$lags[j]))
#' @export
multispectral_h_matrix <- function(CVT, PIT)
  lapply(seq_along(CVT$lags),
         function(j) {
           if (CVT$lags[j]>0) {
             cbind(rep(1,length(PIT)),
                 lagmatrix(CVT$h_closure[[j]](CVT$h_param[[j]])(PIT), CVT$lags[j]))
           } else {
             matrix(1,length(PIT),1)
           }
         })


#' MD regressor matrix
#'
#' Build matrix (or list of matrices) h for the spectral MD test.
#'
#' @param CVT list describing monokernel or multikernel CVT structure
#' @param PIT vector of PIT values (length \eqn{n})
#' @return matrix or list of matrices \eqn{h} (size n x (1+CVT$lags[j]))
#' @export
MD_regressor_matrix <- function(CVT,PIT) {
  switch(CVT$type,
       mono = monospectral_h_matrix(CVT, PIT),
       bi = multispectral_h_matrix(CVT, PIT),
       multi = multispectral_h_matrix(CVT, PIT) )
}

#' Simulation HS
#'
#' Simulate a stylized sample for a hypothetical bank using plain-vanilla historical simulation.
#' Let the unmodelled stochastic volatility be driven by ARMA(1,1) process.
#'
#' @param n length of sample
#' @param cormodel list of \code{ar} and \code{ma} weights, each scalar
#' @return vector of length \code{n}
#' @export
rHS_arma <- function(n,cormodel){
    phi <- cormodel$ar
    psi <- cormodel$ma
    if (!(length(phi)==1) | !(length(psi)==1)) stop ("ARMA(1,1) only")
    if (phi<=0) stop("AR parameter must be positive")
    armadata <- arima.sim(model=cormodel,n=n,n.start=20)
    sigma2 <- (1+2*phi*psi+psi^2)/(1-phi^2)  # true variance
    armadata <- armadata/sqrt(sigma2)  # scale to variance one
    Y <- rbinom(n,1,0.5)
    U <- pnorm(armadata)
    0.5*((1+U)^Y)*(1-U)^(1-Y)
  }

# safeInverse: inverse of Sigma with backup plan if singular
safeInverse <- function(Sigma, cvtname, singular_Sigma_NA)
  tryCatch(solve(Sigma),
           error = function(e) {
                     warning(glue::glue('CVT {cvtname}: Sigma not full rank.\n'))
                     if (singular_Sigma_NA) {
                       NA
                     } else {
                       MASS::ginv(as.matrix(Sigma))
                     }
                   })

# Helper function to enforce single approach to df calculation
dfchisq <- function(invV)
  as.numeric(rankMatrix(invV))

#' Monospectral MD test
#'
#' @param monokernel (list) mono-kernel object (see package help)
#' @param CVT list describing lag structure
#' @param PIT Vector of probability integral transform values.
#' @param h_matrix If provided, then use instead of calling \code{h_closure}.
#' @param informationH (string). Specifies information to incorporate in \eqn{H}.
#' @param singular_Sigma_NA (boolean). If true, then set \eqn{p}-value to \code{NA} when
#'    \eqn{Sigma} is singular.
#' @return \eqn{p}-value of chisq test.
#' @export
monospectral_MDtest <- function(monokernel, CVT, PIT, h_matrix=NULL,
                                   informationH='none', singular_Sigma_NA=FALSE) {
  if (is.null(h_matrix))
      h_matrix <- monospectral_h_matrix(CVT, PIT)
  H <- switch(informationH,
          none = empiricalH(h_matrix),
          partial = partialinfoH(h_matrix),
          full = diag(1+CVT$lags))
  invH <- safeInverse(H, CVT$name, singular_Sigma_NA)
  if (any(is.na(invH))) return(NA)
  dropobs <- is.na(PIT) | is.na(rowSums(h_matrix))
  h <- h_matrix[!dropobs,]
  P <- as.numeric(PIT[!dropobs])
  W <- monokernel$nu(monokernel$support, monokernel$param, standardize=TRUE)(P)
  if (is.degenerateW(W))
    warning(glue::glue('Kernel {monokernel$name}: transformed PIT has degenerate distribution.\n'))
  Y <- as.vector(W %*% h)/length(W)
  stat <- length(W)*mahalanobis(Y, center=0, cov=invH, inverted=TRUE)
  return(1- pchisq(stat, df=dfchisq(invH)))
}

#' Bispectral MD test
#'
#' @param kernlist (list) bi-kernel object (see package help)
#' @param CVT list describing lag structure for each kernel
#' @param PIT Vector of probability integral transform values
#' @param h_list If provided, then use instead of calling \code{h_closure}.
#' @param informationH (string). Specifies information to incorporate in \eqn{H}.
#' @param singular_Sigma_NA (boolean). If true, then set \eqn{p}-value to \code{NA} when
#'    \eqn{Sigma} is singular.
#' @return \eqn{p}-value of chisq test.
#' @export
bispectral_MDtest <- function(bikernel, CVT, PIT, h_list=NULL,
                              informationH='none', singular_Sigma_NA=FALSE) {
  if (is.null(h_list))
      h_list <- multispectral_h_matrix(CVT, PIT)   # works for bispectral too
  H <- switch(informationH,
         none = empiricalH(h_list),
         partial = partialinfoH(h_list),
         full = theoreticalH(CVT))
  rho <- bikernel$correlation(bikernel$support, bikernel$param)
  A_W <- rbind(cbind(matrix(1, nrow=1+CVT$lags[1], ncol=1+CVT$lags[1]),
                     matrix(rho, nrow=1+CVT$lags[1], ncol=1+CVT$lags[2])),
               cbind(matrix(rho, nrow=1+CVT$lags[2], ncol=1+CVT$lags[1]),
                     matrix(1, nrow=1+CVT$lags[2], ncol=1+CVT$lags[2])))
  invSigma <- safeInverse(A_W * H, CVT$name, singular_Sigma_NA)
  if (any(is.na(invSigma))) return(NA)
  keepobs <- !is.na(PIT) & !is.na(rowSums(do.call(cbind,h_list)))
  P <- as.numeric(PIT[keepobs])
  W_list <- apply_nu_to_PIT(bikernel, P, standardize = TRUE)
  if (any(is.degenerateW(W_list)))
    warning(glue::glue('Kernel {bikernel$name}: transformed PIT has degenerate distribution.\n'))

  Y <- as.vector(c(W_list[[1]] %*% h_list[[1]][keepobs,],
                   W_list[[2]] %*% h_list[[2]][keepobs,]))/length(P)
  stat <- length(P)*mahalanobis(Y, center=0, cov=invSigma, inverted=TRUE)
  return(1- pchisq(stat, df=dfchisq(invSigma)))
}

#' Multispectral SUR MD test
#'
#' @param kernlist (list) multi- or bi-kernel object (see package help)
#' @param CVT list describing lag structure for each kernel
#' @param PIT Vector of probability integral transform values
#' @param h_list If provided, then use instead of calling \code{h_closure}.
#' @param GW (boolean) If true, then run Giacomini-White (OLS) in place of SUR.
#' @param singular_Sigma_NA (boolean). If true, then set \eqn{p}-value to \code{NA} when
#'    \eqn{Sigma} is singular.
#' @return \eqn{p}-value of chisq test.
#' @export
multispectral_SUR_MDtest <- function(kernlist, CVT, PIT, h_list=NULL, GW=FALSE,
                                     singular_Sigma_NA=FALSE, npval=0) {
  stopifnot(length(CVT$lags)>=2)
  if (is.null(h_list))
    h_list <- multispectral_h_matrix(CVT, PIT)   # works for bispectral too
  keepobs <- !is.na(PIT) & !is.na(rowSums(do.call(cbind,h_list)))
  P <- as.numeric(PIT[keepobs])
  # Form block diagonal matrix of regressors
  X <- map(h_list, ~.x[keepobs,]) %>% .bdiag
  # Independent variable first in list form
  W_list <- apply_nu_to_PIT(kernlist, P, standardize = TRUE) %>% map(as.matrix)
  if (any(is.degenerateW(W_list)))
    warning(glue::glue('Kernel {kernlist$name}: transformed PIT has degenerate distribution.\n'))
  Y <- do.call(rbind, W_list)
  rhoMat <- multispectral_correlation(kernlist)
  if (GW) {
    Omega <- kronecker(rhoMat, Diagonal(length(P)))
    invSigma <- safeInverse(Matrix::t(X) %*% Omega %*% X, CVT$name, singular_Sigma_NA)
    Xt <- Matrix::t(X)
  } else {
    invOmega <- kronecker(solve(rhoMat), Diagonal(length(P)))
    invSigma <- safeInverse(Matrix::t(X) %*% invOmega %*% X, CVT$name, singular_Sigma_NA)
    Xt <- Matrix::t(X) %*% invOmega
  }
  if (any(is.na(invSigma))) return(NA)
  XtY <- Xt %*% Y
  stat <- as.numeric(Matrix::t(XtY) %*% invSigma %*% XtY)
  if (npval==0) {
     return(1- pchisq(stat, df=dfchisq(invSigma)))
  } else {
     rpval <- function() {
       W_list <- apply_nu_to_PIT(kernlist, runif(length(P)), standardize = TRUE) %>% map(as.matrix)
       XtY <- Xt %*% do.call(rbind, W_list)
       as.numeric(Matrix::t(XtY) %*% invSigma %*% XtY)
     }
     Femp <- replicate(npval, rpval()) %>% ecdf()
     return(1-Femp(stat))
  }
}

#' Spectral MD test
#'
#' @param kernel (list) kernel object (see package help)
#' @param CVT (list) CVT object (see package help)
#' @param PIT Vector of probability integral transform values
#' @param h_list If provided, then use instead of calling \code{h_closure}.
#' @param informationH (string). Specifies information to incorporate in \eqn{H}.
#' @param singular_Sigma_NA (boolean). If true, then set \eqn{p}-value to \code{NA} when
#'    \eqn{Sigma} is singular.
#' @return \eqn{p}-value of chisq test.
#' @export
spectral_MDtest <- function(kernel, CVT, PIT, h_list=NULL,
                            informationH='none', singular_Sigma_NA=FALSE) {
  stopifnot(kernel$type==CVT$type)
  switch(kernel$type,
         mono = monospectral_MDtest(kernel, CVT, PIT, h_list, informationH, singular_Sigma_NA),
         bi = bispectral_MDtest(kernel, CVT, PIT, h_list, informationH, singular_Sigma_NA),
         multi = stop("Multi-spectral MD test not yet implemented."))
}


