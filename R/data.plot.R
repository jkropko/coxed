#' Plot the histograms of simulated data
#'
#' This function is called by \code{\link[coxed]{survsim.plot}} and is not intended to be used by itself.
#' @param data A data frame of simulated duration data in which the duration variable is called y.
#' Generally, this data frame is taken from the \code{data} attribute of the \code{\link[coxed]{sim.survdata}}
#' function
#' @param xb A vector of the linear predictors, combining the covariates and coefficients from the simulation.
#' Generally, this data frame is taken from the \code{xb} attribute of the \code{\link[coxed]{sim.survdata}}
#' function
#' @param bins The number of bins to draw in each histogram
#' @details This function produces three histograms, one of the simulated durations, one of the linear predictors,
#' and one of the exponentiated linear predictors. It uses \code{\link[ggplot2]{facet_wrap}} within
#' a call to \code{\link[ggplot2]{ggplot}} to arrange the plots.
#' @return A figure of class \code{"gg"} and \code{"ggplot"}
#' @export
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @seealso \code{\link[coxed]{survsim.plot}}, \code{\link[coxed]{sim.survdata}}
#' @examples
#' simdata <- sim.survdata(N=1000, T=100, num.data.frames=1)
#' data.plot(simdata$data, simdata$xb)
data.plot <- function(data, xb, bins=30){
  d <- data.frame(quantity = "Simulated durations", value = data$y)
  d <- rbind(d, data.frame(quantity = "Linear predictor",
                           value = xb))
  d <- rbind(d, data.frame(quantity = "Exponentiated linear predictor",
                           value = exp(xb)))
  g <- ggplot2::ggplot(d, ggplot2::aes(x = value)) +
       ggplot2::geom_histogram(bins=bins) +
       ggplot2::facet_wrap(~ quantity, scales = "free") +
       ggplot2::xlab("Value") +
       ggplot2::ylab("Frequency") +
       ggplot2::ggtitle("Histograms of Simulated Data")
  return(g)
}
