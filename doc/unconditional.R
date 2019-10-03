## ----echo=FALSE, message=FALSE-------------------------------------------
library(spectralBacktest)
library(purrr)
library(dplyr)
library(knitr)
set.seed(543)


## ----fig.width=6, fig.height=4-------------------------------------------

n <- 750
alfbet <- 0.8

PIT1 <- runif(n)
PIT2 <- rbeta(n,alfbet,alfbet)

hist(PIT2)


## ----warning=FALSE-------------------------------------------------------

alpha1 <- 0.95
alpha_star <- 0.99
alpha2 <- 0.995

BIN <- list( name = 'Binomial score at 99%',
             type = 'mono',
             nu = nu_discrete,
             support = alpha_star,
             param = 1 )

ZU3 <- list( name = 'Discrete Uniform 3',
             type = 'mono',
             nu = nu_discrete,
             support = c(alpha1, alpha_star, alpha2),
             param = c(1, 1, 1) )

ZE <- list( name = 'Epanechnikov',
            type = 'mono',
            nu = nu_epanechnikov,
            support = c(alpha1, alpha2),
            param = NULL )

ZU <- list( name = 'Uniform',
            type = 'mono',
            nu = nu_uniform,
            support = c(alpha1, alpha2),
            param = NULL )

ZLp <- list( name = 'LinearUp',
            type = 'mono',
            nu = nu_linear,
            support = c(alpha1, alpha2),
            param = 1 )
ZLn <- list( name = 'LinearDown',
            type = 'mono',
            nu = nu_linear,
            support = c(alpha1, alpha2),
            param = -1 )


ZLL <- list( name = 'linear/Linear',
             type = 'bi',
             nu = list(nu_linear, nu_linear),
             correlation = rho_linear_linear,
             support = c(alpha1, alpha2),
             param = list(-1,1) )

ZAE <- list( name = 'Arcsin/Epanechnikov',
             type = 'bi',
             nu = list(nu_arcsin, nu_epanechnikov),
             correlation = rho_arcsin_epanechnikov,
             support = c(alpha1, alpha2),
             param = list(NULL, NULL) )

PNS <- list( name = 'Probitnormal score',
             type = 'bi',
             nu = list(nu_probitnormal, nu_probitnormal),
             correlation = rho_probitnormal,
             support = c(alpha1, alpha2),
             param = list(1, 2) )

Pearson3 <- list(name = 'Pearson',
                 type = 'multi',
                 nu = nu_pearson,
                 correlation = rho_pearson,
                 support=NULL,
                 param=list(alpha1, alpha_star, alpha2))

# gather the tests into a list and execute!
kernlist <- list(BIN=BIN, Pearson3=Pearson3, ZU3=ZU3,
                 ZU=ZU, ZE=ZE, ZLp=ZLp, ZLn=ZLn,
                 ZLL=ZLL, ZAE=ZAE, PNS=PNS)

pval <- map_df(list(PIT1, PIT2), 
               function(P) lapply(kernlist, function(kern) spectral_Ztest(kern,P))) %>%
        mutate(PIT=c('Uniform',sprintf('Beta(%0.2f,%0.2f)',alfbet,alfbet))) %>%
        select('PIT', everything())

kable(pval, digits=4)


## ----warning=FALSE-------------------------------------------------------

n <- 500
a <- b <- 1  # if a==b==1, we are assessing size of test
Npf <- 5000  # number of portfolios
rpval <- function(kernellist) {
  P <- rbeta(n,a,b)
  map(kernellist, ~spectral_Ztest(.x,P))
}

df <- rerun(Npf, rpval(kernlist)) %>%  map_df(`[`,names(kernlist))
rejectrate <- summarize_all(df, list(function(x) mean(x<=0.05)))
kable(rejectrate, digits=4, caption='Size. Frequency of test rejections at 5% level')


## ----warning=FALSE-------------------------------------------------------

a <- b <- alfbet  
df <- rerun(Npf, rpval(kernlist)) %>%  map_df(`[`,names(kernlist))
rejectrate <- summarize_all(df, list(function(x) mean(x<=0.05)))
kable(rejectrate, digits=4, caption='Power. Frequency of test rejections at 5% level')


