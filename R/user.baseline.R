#' Calculating baseline functions from a user-specified baseline hazard function
#'
#' This function is called by \code{\link[coxed]{sim.survdata}} and is not intended to be used by itself.
#' @param user.fun A user-specified R function with one argument, representing time, that outputs the baseline hazard function
#' @param T The latest time point during which an observation may fail. Failures can occur
#' as early as 1 and as late as T
#' @details \code{user.baseline} takes a function as a user-specified baseline hazard which must have
#' only one argument: time.  \code{user.baseline} approximates the cumulative baseline hazard by taking the
#' cumulative sum of the user-specified hazard function. It calculates the survivor function by exponentiating
#' the cumulative baseline hazard time -1, the baseline failure CDF by subtracting the survivor function from 1,
#' and it approximates the baseline failure PDF by taking the first difference of the failure CDF.
#' survivor function, and baseline failure-time PDF and CDF.
#' @return A data frame with five columns representing time from 1 to T, and the user-specified baseline hazard,
#' cumulative hazard, survivor function, failure PDF and failure CDF at each time point.
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @export
#' @examples
#' ## Writing the hazard to be lognormal with mean of 50, sd of 10
#' my.hazard <- function(t){
#'      dnorm((log(t) - log(50))/log(10)) /
#'           (log(10)*t*(1 - pnorm((log(t) - log(50))/log(10))))
#' }
#' lognormal.functions <- user.baseline(my.hazard, 100)
#' summary(lognormal.functions)
#'
#' #A customized user-specified hazard
#' sine.squared.hazard <- user.baseline(function(t) sin(t/25)^2, 30)
#' summary(sine.squared.hazard)
user.baseline <- function(user.fun, T){
     baseline <- data.frame(time=1:T)
     baseline <- dplyr::mutate(baseline,
                               hazard = user.fun(time),
                               cum.hazard = cumsum(hazard),
                               survivor = exp(-cum.hazard),
                               failure.CDF = 1 - survivor,
                               failure.PDF = c(0, diff(failure.CDF)))
     baseline <- dplyr::select(baseline, -cum.hazard)
     return(baseline)
}
