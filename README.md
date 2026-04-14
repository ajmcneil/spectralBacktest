# spectralBacktest
This `R` package implements spectral backtests of forecast distributions. The methodology is described in two papers by Michael B. Gordy and Alexander J. McNeil:

* *Spectral backtests of forecast distributions with application to risk management*, published in *Journal of Banking and Finance* **116**, July 2020.  [JBF](https://doi.org/10.1016/j.jbankfin.2020.105817) 
* *Spectral backtests unbounded and folded*, conditionally accepted to *Journal of Financial Econometrics*. [Preprint on FEDS](https://doi.org/10.17016/FEDS.2024.060)

To install:
> devtools::install_github("ajmcneil/spectralBacktest")

To install with vignettes (SLOWER!):
> devtools::install_github("ajmcneil/spectralBacktest", build_vignettes = TRUE)
