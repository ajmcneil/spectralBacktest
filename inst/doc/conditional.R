## ----echo=FALSE, message=FALSE-------------------------------------------
library(spectralBacktest)
library(purrr)
library(knitr)
set.seed(543)


## ----fig.width=7, fig.height=5, echo=TRUE, results='asis'----------------

n <- 500

P <- rHS_arma(n, list(ar=0.35, ma=-0.1))
acf(P,lag=5)
acf(abs(2*P-1),lag=5)


## ----warning=FALSE-------------------------------------------------------

alpha1 <- 0.95
alpha_star <- 0.99
alpha2 <- 0.995


ZU <- list( name = 'Uniform',
            type = 'mono',
            nu = nu_uniform,
            support = c(alpha1, alpha2),
            param = NULL )

ZLL <- list( name = 'Linear/Linear',
             type = 'bi',
             nu = list(nu_linear, nu_linear),
             correlation = rho_linear_linear,
             support = c(alpha1, alpha2),
             param = list(-1, 1) )

PNS <- list( name = 'Probitnormal score',
             type = 'bi',
             nu = list(nu_probitnormal, nu_probitnormal),
             correlation = rho_probitnormal,
             support = c(alpha1, alpha2),
             param = list(1L, 2L) )

# gather the tests into a list and execute!
kernlist <- list(ZU=ZU, ZLL=ZLL, PNS=PNS)
pval_Z <- map_df(kernlist, function(kern) spectral_Ztest(kern,P))
kable(pval_Z, digits=5, caption='p-values of tests of unconditional coverage')



## ----warning=TRUE--------------------------------------------------------

# Monospectral
CVT <- list( name = 'Power 4',
             type = 'mono',
             h_closure = CVT_PtildePower,
             h_param = 4,
             lags = 4L )

# Bispectral.  Note that each element of lags is a nonnegative integer.
CVT2 <- list( name = 'Power4/Power 0.5',
              type = 'bi',
              h_closure = list(CVT_PtildePower, CVT_PtildePower),
              correlation = rho_PtildePower,
              h_param = list(4,1/2),
              lags = c(4L,2L) )

MDlist <- list(CVT, CVT2, CVT2)  # same length as kernlist, defined earlier

infoinH <- list('none', 'partial', 'full')
pval_MD <- map_df(infoinH, function(s) 
                       q <- c(Information=s,
                              map2(kernlist, MDlist, 
                              ~spectral_MDtest(.x, .y, P, informationH=s))))

kable(pval_MD, digits=5, caption='p-values of MD tests of conditional coverage')




