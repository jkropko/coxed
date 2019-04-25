#' Marginal changes in expected duration
#'
#' This function reports summary statistics for the changes in duration that occur when
#' comparing observations in \code{newdata2} to observations in \code{newdata} as specified
#' in the original call to \code{\link[coxed]{coxed}}.
#' @param object The output from \code{\link[coxed]{coxed}}. If neither \code{newdata} or
#' \code{newdata2} are \code{NULL}, \code{coxed} is being calculate a marginal
#' effect and the class of the \code{coxed} output will be "\code{coxedMargin}", so
#' this function will be called by the generic \code{summary} function
#' @param stat Either \code{"mean"} or \code{"median"}
#' @param ... For future methods
#' @details \code{coxed} calculates a vector of expected durations for every observation
#' in \code{newdata2} and \code{newdata} (these data frames must have the same number of
#' observations, and corresponding rows are assumed to refer to the same observation), and
#' subtracts the expected duration from \code{newdata2} from the expected duration from
#' \code{newdata}.  That generates a vector of marginal changes in duration for every
#' observation.  To generalize a finding across observations, either the mean (if \code{type="mean"})
#' or the median (if \code{type="median"}) of these differences is reported, along with the
#' mean/median of the expected durations from each of the two covariate profiles.
#' If \code{bootstrap=TRUE} in the call to \code{coxed} then a bootstrapped standard error
#' and confidence interval is reported for each of these quantities as well.
#' @return A data.frame containing the mean or median duration for each of \code{newdata2} and
#' \code{newdata} and the differences in these durations across the two covariate profiles.
#' If \code{bootstrap=TRUE} in the call to \code{coxed}, the data frame also includes
#' the bootstrapped standard error and confidence interval for these quantities.
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @seealso \code{\link[coxed]{coxed}}
#' @export
#' @examples
#' mv.surv <- Surv(martinvanberg$formdur, event = rep(1, nrow(martinvanberg)))
#' mv.cox <- coxph(mv.surv ~ postel + prevdef + cont + ident + rgovm + pgovno + tpgovno +
#'      minority, method = "breslow", data = martinvanberg)
#' summary(mv.cox)
#'
#' me <- coxed(mv.cox, method="npsf", bootstrap = FALSE,
#'             newdata = dplyr::mutate(martinvanberg, rgovm=0),
#'             newdata2 = dplyr::mutate(martinvanberg, rgovm=1.24))
#' summary(me, stat="mean")
#' summary(me, stat="median")
summary.coxedMargin <- function(object, stat="mean", ...) {
     stopifnot(inherits(object, "coxedMargin"))
     if(!(stat %in% c("mean", "median"))) stop("stat must be one of mean or median")
     if(stat=="mean"){
          if(is.null(object$mean1)){
               r <- colMeans(object$exp.dur, na.rm=TRUE)
               names(r) <- c("newdata2", "newdata", "difference")
          } else{
               r <- rbind(object$mean2, object$mean1, object$mean.diff)
               rownames(r) <- c("newdata2", "newdata", "difference")
          }
     }
     if(stat=="median"){
          if(is.null(object$mean1)){
               r <- apply(object$exp.dur, 2, FUN=function(x){
                    median(x, na.rm=TRUE)
               })
               names(r) <- c("newdata2", "newdata", "difference")
          } else {
               r <- rbind(object$median2, object$median1, object$median.diff)
               rownames(r) <- c("newdata2", "newdata", "difference")
          }
     }
     r <- round(r, 3)
     return(r)
}
