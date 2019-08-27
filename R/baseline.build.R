#' Generate simulated baseline hazard, cumulative hazard, survival, failure PDF, and failure CDF functions
#'
#' This function is called by \code{\link[coxed]{sim.survdata}} and is not intended to be used by itself.
#' @param T The latest time point during which an observation may fail. Failures can occur
#' as early as 1 and as late as T
#' @param knots The number of points to draw while using the flexible-hazard method to generate hazard functions (default is 8).
#' Ignored if \code{hazard.fun} is not \code{NULL}.
#' @param spline If \code{TRUE} (the default), a spline is employed to smooth the generated cumulative baseline hazard, and if \code{FALSE}
#' the cumulative baseline hazard is specified as a step function with steps at the knots. Ignored if \code{hazard.fun} is not \code{NULL}
#' @details This function employs the flexible hazard method described in Harden and Kropko (2018) to generate a baseline
#' failure CDF: it plots points at (0, 0) and (\code{T}+1, 1), and it plots \code{knots} additional points with x-coordinates drawn uniformly
#' from integers in [2, \code{T}] and y-coordinates drawn from U[0, 1]. It sorts these coordinates in ascending order
#' (because a CDF must be non-decreasing) and if \code{spline=TRUE} it fits a spline using Hyman’s (1983) cubic smoothing function to preserve the CDF’s monotonicity.
#' Next it constructs the failure-time PDF by computing the first differences of the CDF at each time point.
#' It generates the survivor function by subtracting the failure CDF from 1. Finally, it computes the baseline hazard by dividing the failure PDF by the survivor function.
#' @return A data frame containing every potential failure time and the baseline failure PDF,
#' baseline failure CDF, baseline survivor function, and baseline hazard function at each time point.
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @export
#' @references Harden, J. J. and Kropko, J. (2018). Simulating Duration Data for the Cox Model.
#' \emph{Political Science Research and Methods} \url{https://doi.org/10.1017/psrm.2018.19}
#'
#' Hyman, J. M. (1983) Accurate monotonicity preserving cubic interpolation. \emph{SIAM J. Sci. Stat. Comput.} \strong{4}, 645–654.
#' @seealso \code{\link[stats]{splinefun}}, \code{\link[coxed]{sim.survdata}}
#' @examples
#' baseline.functions <- baseline.build(T=100, knots=8, spline=TRUE)
baseline.build <- function(T=100, knots = 8, spline = TRUE){
  time <- 1:(T+1)
  k <- c(1,sort(sample(time[2:T], size=knots, replace=FALSE)), (T+1))
  heights <- c(0, sort(runif(knots)), 1)
  tk <- merge(data.frame(time), data.frame(time=k, heights),
              by="time", all = !spline)
  MonotonicSpline <- stats::splinefun(x = tk$time, y = tk$heights,
                               method = "hyman")
  bl.failure.CDF <- MonotonicSpline(time)
  baseline <- data.frame(time = time[-(T+1)],
                         failure.PDF = diff(bl.failure.CDF),
                         failure.CDF = bl.failure.CDF[-1],
                         survivor = abs(1 - bl.failure.CDF[-1]))
  baseline$hazard <-baseline$failure.PDF/(1 - bl.failure.CDF[-(T+1)])
  return(baseline)
}
