<!-- README.md is generated from README.Rmd. Please edit that file -->
Duration-Based Quantities of Interest and Simulation Methods for the Cox Proportional Hazards Model
===================================================================================================

The Cox proportional hazards model (implemented in R with `coxph()` in the `survival` package or with `cph()` in the `rms` package) is one of the most frequently used estimators in duration (survival) analysis. Because it is estimated using only the observed durations' rank ordering, typical quantities of interest used to communicate results of the Cox model come from the hazard function (e.g., hazard ratios or percentage changes in the hazard rate). These quantities are substantively vague and difficult for many audiences of research to understand.

The `coxed` package introduces a suite of methods to address these problems. The package allows researchers to calculate duration-based quantities from Cox model results, such as the expected duration (or survival time) given covariate values and marginal changes in duration for a specified change in a covariate. These duration-based quantities often match better with researchers' substantive interests and are easily understood by most readers.

In addition, no standard method exists for simulating durations directly from the Cox model's data generating process because it does not assume a distributional form for the baseline hazard function. The `coxed` package also contains functions to simulate general duration data that does not rely on an assumption of any particular parametric hazard function.

Examples
--------

Please see our two vignettes: `coxed.Rmd` and `simulating_survival_data.Rmd` for detailed instructions with examples about how to use the `coxed()` and `sim.survdata()` functions.

``` r
library(coxed)
```

### A simple example of the `coxed()` function

First we replicate the Cox model from Martin and Vanberg (2003):

``` r
mv.surv <- Surv(martinvanberg$formdur, event = rep(1, nrow(martinvanberg)))
mv.cox <- coxph(mv.surv ~ postel + prevdef + cont + ident + rgovm + pgovno + 
                     tpgovno + minority, method = "breslow", data = martinvanberg)
summary(mv.cox)
#> Call:
#> coxph(formula = mv.surv ~ postel + prevdef + cont + ident + rgovm + 
#>     pgovno + tpgovno + minority, data = martinvanberg, method = "breslow")
#> 
#>   n= 203, number of events= 203 
#> 
#>              coef exp(coef) se(coef)       z Pr(>|z|)    
#> postel   -0.57665   0.56177  0.16862  -3.420 0.000627 ***
#> prevdef  -0.10000   0.90484  0.22987  -0.435 0.663543    
#> cont      1.10047   3.00556  0.23970   4.591 4.41e-06 ***
#> ident     0.14579   1.15695  0.11859   1.229 0.218938    
#> rgovm    -0.21312   0.80806  0.12036  -1.771 0.076625 .  
#> pgovno    1.19057   3.28894  0.12405   9.598  < 2e-16 ***
#> tpgovno  -0.43189   0.64928  0.03476 -12.425  < 2e-16 ***
#> minority -0.42759   0.65208  0.20793  -2.056 0.039740 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>          exp(coef) exp(-coef) lower .95 upper .95
#> postel      0.5618     1.7801    0.4037    0.7818
#> prevdef     0.9048     1.1052    0.5766    1.4198
#> cont        3.0056     0.3327    1.8789    4.8079
#> ident       1.1570     0.8643    0.9170    1.4597
#> rgovm       0.8081     1.2375    0.6382    1.0231
#> pgovno      3.2889     0.3040    2.5791    4.1942
#> tpgovno     0.6493     1.5402    0.6065    0.6951
#> minority    0.6521     1.5336    0.4338    0.9801
#> 
#> Concordance= 0.903  (se = 0.025 )
#> Rsquare= 0.745   (max possible= 1 )
#> Likelihood ratio test= 277.2  on 8 df,   p=0
#> Wald test            = 218.1  on 8 df,   p=0
#> Score (logrank) test = 279.3  on 8 df,   p=0
```

To see predicted durations from the Cox model, place the Cox model output as the first argument of `coxed()`:

``` r
ed1 <- coxed(mv.cox, method="npsf")
```

There are a number of uses of the `coxed()` output. First, the predicted durations for each individual observation are stored in the `exp.dur` attribute:

``` r
head(ed1$exp.dur)
#>     exp.dur
#> 1 34.620664
#> 2 31.639975
#> 3 39.509093
#> 4 14.717544
#> 5  2.473799
#> 6 47.347670
```

The `summary()` function, when applied to `coxed`, reports either the mean or median duration in the estimation sample, depending on the option specified with `stat`:

``` r
summary(ed1, stat="mean")
#>  mean 
#> 25.18
summary(ed1, stat="median")
#> median 
#> 19.121
```

The predicted mean duration of government negotiations is 25.18 days, and the predicted median duration is 19.12 days.

In addition to reporting the mean and median duration, the NPSF version of `coxed()` provides estimates of the cumulative baseline hazard function and the baseline survivor function in the data. These functions are stored as a data frame in the `baseline.functions` attribute.

``` r
head(ed1$baseline.functions)
#>   time        cbh  survivor
#> 1    1 0.01631230 0.9838200
#> 2    2 0.03587341 0.9647624
#> 3    3 0.04746259 0.9536462
#> 4    4 0.07651484 0.9263392
#> 5    5 0.09169880 0.9123799
#> 6    6 0.12573075 0.8818522
```

### Simulating a single dataset

To generate a survival dataset, use the `sim.survdata()` function, with the `num.data.frames` argument set to 1. Here we generate a single survival dataset with 1000 observations, in which durations can fall on any integer between 1 and 100:

``` r
simdata <- sim.survdata(N=1000, T=100, num.data.frames=1)
```
