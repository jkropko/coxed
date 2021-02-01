#' Calculating a simulated marginal effect
#'
#' This function is called by \code{\link[coxed]{sim.survdata}} and is not intended to be used by itself.
#' @param baseline The baseline hazard functions, output by \code{\link[coxed]{baseline.build}}
#' @param xb The simulated data, output by \code{\link[coxed]{generate.lm}}
#' @param covariate Specification of the column number of the covariate in the \code{X} matrix for which to generate a simulated marginal effect (default is 1).
#' The marginal effect is the difference in expected duration when the covariate is fixed at a high value and the expected duration when the covariate is fixed
#' at a low value
#' @param low The low value of the covariate for which to calculate a marginal effect
#' @param high The high value of the covariate for which to calculate a marginal effect
#' @param compare The statistic to employ when examining the two new vectors of expected durations (see details for \code{\link[coxed]{sim.survdata}}).  The default is \code{median}
#' @details The idea is to simulate a marginal change in duration so that researchers can compare the performance of
#' estimators of this statistic using simulated data.
#'
#' The function calculates simulated durations for each observation conditional on a baseline hazard function
#' and exogenous covariates and coefficients.  The \code{covariate} argument specifies the variable in the X matrix to
#' vary so as to measure the marginal effect.  First the covariate is set to the value specified in \code{low} for all
#' observations, then to the value specified in \code{high} for all observations.  Given each value, new durations are
#' drawn. The durations when the covariate equals the low value are subtracted from the durations when the covariate
#' equals the high value.  The marginal effect is calculated by employing the statistic given by \code{compare}, which
#' is \code{median} by default.
#' @return A list with three items:
#' \tabular{ll}{
#' \code{marg.effect} \tab A scalar containing the simulated marginal effect\cr
#' \code{data.low} \tab The durations and covariates when the covariate of interest is set to the low value \cr
#' \code{data.high} \tab The durations and covariates when the covariate of interest is set to the high value \cr
#' }
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @seealso \code{\link[coxed]{baseline.build}}, \code{\link[coxed]{generate.lm}}, \code{\link[coxed]{sim.survdata}}
#' @export
#' @examples
#' T <- 100
#' N <- 1000
#' X <- as.matrix(data.frame(X1=rnorm(N), X2=rnorm(N), X3=rnorm(N)))
#' beta <- as.matrix(rnorm(3))
#' baseline <- baseline.build(T=T, knots=8, spline=TRUE)
#' xb <- generate.lm(baseline, X=X, beta=beta, N=N, censor=.1, type="none")
#' me <- make.margeffect(baseline, xb, covariate=1, low=0, high=1)
#' me$marg.effect
make.margeffect <- function(baseline, xb, covariate=1, low=0, high=1, compare=median){

        if(xb$tvc){
                X0 <- dplyr::select(xb$data, -id, -failed, -start, -end)
                X1 <- dplyr::select(xb$data, -id, -failed, -start, -end)
        } else {
                X0 <- dplyr::select(xb$data, -y)
                X1 <- dplyr::select(xb$data, -y)
        }

        X0[,covariate] <- low
        X1[,covariate] <- high

        beta <- xb$beta

        if(nrow(as.matrix(beta)) == 1){
                XB0 <- as.matrix(X0)%*% beta
                survival <- t(sapply(XB0, FUN=function(x){baseline$survivor^exp(x)}, simplify=TRUE))
        } else {
                XB0 <- apply(beta, 1, FUN = function(b){
                        return(as.matrix(X0) %*% b)
                })
                survival <- t(apply(XB0, 1, FUN=function(x){baseline$survivor^exp(x)}))
        }
        y0 <- apply(survival, 1, FUN=function(x){
                which.max(diff(x < runif(1)))
        })

        if(nrow(as.matrix(beta)) == 1){
                XB1 <- as.matrix(X1)%*%beta
                survival <- t(sapply(XB1, FUN=function(x){baseline$survivor^exp(x)}, simplify=TRUE))
        } else {
                XB1 <- apply(beta, 1, FUN = function(b){
                        return(as.matrix(X1) %*% b)
                })
                survival <- t(apply(XB1, 1, FUN=function(x){baseline$survivor^exp(x)}))
        }
        y1 <- apply(survival, 1, FUN=function(x){
                which.max(diff(x < runif(1)))
        })

        marg.effect <- compare(y1 - y0)
        data.low <- list(x = X0, y = y0)
        data.high <- list(x = X1, y = y1)

        return(list(marg.effect = marg.effect,
                    data.low = data.low,
                    data.high = data.high))
}
