#' Generate right-censoring to be dependent on the data
#'
#' This function is called by \code{\link[coxed]{sim.survdata}} and is not intended to be used by itself.
#' @param x A matrix or data frame containing covariates
#' @param censor The proportion of observations to designate as being right-censored
#' @details The purpose of this function is to efficiently generate indicators for whether or not an observation in simulated duration
#' data is right-censored. In this case, whether or not an observation is right-censored depends on the covariates in the data.
#'
#' This function randomly draws new coefficients, one for every column of \code{x}, from the normal distribution with mean 0 and
#' standard deviation of 0.1.  It uses these new coefficients to build a linear predictor, to which is added a disturbance term which
#' is also drawn from N(0,.1). A fixed proportion, given by \code{censor}, of the observations with the highest values of this linear
#' predictor are set to be \code{TRUE} and the others are set to \code{FALSE}.
#' @return A vector of logical values
#' @examples Xdata <- matrix(rnorm(300), 100, 3)
#' censor.x(Xdata, .1)
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @seealso \code{\link[coxed]{sim.survdata}}
#' @export
censor.x <- function(x, censor=.1){
  beta.cen <- as.matrix(rnorm(ncol(x), mean=0, sd=.1))
  xb.cen <- as.matrix(x)%*%beta.cen + rnorm(nrow(x), mean=0, sd=.1)
  cen <- (xb.cen > quantile(xb.cen, (1-censor)))
  return(cen)
}
