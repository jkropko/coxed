#' Plot the simulated baseline functions
#'
#' This function is called by \code{\link[coxed]{survsim.plot}} and is not intended to be used by itself.
#' @param baseline A data frame containing five variables: time, and the values of the baseline failure PDF,
#' baseline failure CDF, baseline survivor function, and baseline hazard function at each time point.
#' Generally, this data frame is taken from the \code{baseline} attribute of the \code{\link[coxed]{sim.survdata}}
#' function
#' @details This function reshapes the data for easy faceting with \code{\link[ggplot2]{facet_wrap}} within
#' a call to \code{\link[ggplot2]{ggplot}}. Each function is plotted on the y-axis and time is plotted on
#' the x-axis using \code{\link[ggplot2]{geom_line}}
#' @return A figure of class \code{"gg"} and \code{"ggplot"}
#' @export
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @seealso \code{\link[coxed]{survsim.plot}}, \code{\link[coxed]{sim.survdata}}
#' @examples
#' simdata <- sim.survdata(N=1000, T=100, num.data.frames=1)
#' baseline.plot(simdata$baseline)
baseline.plot <- function(baseline){
     baseline <- tidyr::gather(baseline, failure.PDF, failure.CDF,
                               survivor, hazard, key="type", value="value")
     baseline$type <- factor(baseline$type,
                             levels = c("failure.PDF", "failure.CDF", "survivor", "hazard"),
                             labels = c("Failure PDF", "Failure CDF", "Survivor", "Hazard"))
     g <- ggplot2::ggplot(baseline, ggplot2::aes(x=time, y=value)) +
          ggplot2::geom_line() +
          ggplot2::facet_wrap(~ type, scales="free") +
          ggplot2::xlab("Time") +
          ggplot2::ylab("Survival Model Function") +
          ggplot2::ggtitle("Simulated Baseline Functions")
     return(g)
}
