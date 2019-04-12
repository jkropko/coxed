#' Generate predicted ranks for new observations given a new covariate profile
#'
#' This function is called by \code{\link[coxed]{coxed}} when \code{method="gam"} and new data
#' are specified, and is not intended to be used by itself.
#' @param x A vector of linear predictors for the estimation sample
#' @param v A vector of linear predictors for the new data
#' @param ties.method A character string specifying how ties are treated,
#' see ‘Details’ in the documentation for \code{\link[base]{rank}}; can be
#' abbreviated.
#' @param warn If \code{TRUE} the function warns the user when linear predictors in the new data are greater than
#' or less than all of the linear predictors in the estimation sample
#' @details The purpose of \code{rank.predict} is to determine for a single new observation what the rank of that
#' observation's linear predictor would have been had the observation been in the original estimation sample.
#' It calculates the predicted rank by appending the new observation to the vector of linear predictors for the
#' estimation sample and calculating the \code{\link[base]{rank}} of the observation in the new vector. If
#' the new data contain more than one observation, \code{rank.predict} calculates the predicted rank for each
#' observation independently, without taking the other observations in the new data into account.
#'
#' Any observation with a linear predictor less than the minimum linear predictor in the estimation sample will have a predicted
#' rank of 1; any observation with a linear predictor greater than the maximum linear predictor in the estimation
#' sample will have a predicted rank of \code{length(v)}. If either condition exists, the function provides a warning.
#' @return A numeric vector containing the predicted ranks for the observations in \code{x}.
#' @seealso \code{\link[coxed]{coxed}}, \code{\link[base]{rank}}
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @export
#' @examples
#' estimationLPs <- rnorm(20)
#' cbind(estimationLPs, rank(estimationLPs))
#' newLPs <- rnorm(5)
#' newLP.rank <- rank.predict(x=newLPs, v=estimationLPs)
#' cbind(newLPs, newLP.rank)
rank.predict <- function(x, v, ties.method="random", warn=TRUE){
     r <- sapply(x, FUN=function(t){
          v2 <- rbind(data.frame(v=v, kind=0), data.frame(v=t, kind=1))
          v2$v.rank <- rank(v2$v, ties.method=ties.method)
          v2 <- v2[v2$kind==1,]
          return(v2$v.rank)
     })
     c1 <- sum(r==1)
     c2 <- sum(r > length(v))
     paste1 <- paste(c("New data contain", c1,
                       "observations with linear predictors less than all linear predictors in the estimation sample. These observations will all have the same predicted duration"),
                     collapse=" ")
     paste2 <- paste(c("New data contain", c2,
                       "observations with linear predictors greater than all linear predictors in the estimation sample. These observations will all have the same predicted duration"),
                     collapse=" ")
     if(warn & min(r)==1) warning(paste1)
     if(warn & max(r)>length(v)) warning(paste2)
     return(r)
}
