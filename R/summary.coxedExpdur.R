#' The mean or median expected duration
#'
#' This function takes the output of \code{\link[coxed]{coxed}} and calculates the
#' mean or median of the expected durations across the observations for which durations
#' are estimated.
#' @param object The output from \code{\link[coxed]{coxed}}. If \code{newdata2=NULL}, so that
#' \code{coxed} is being used to predict duration instead of to calculate a marginal
#' effect, then the class of the \code{coxed} output will be "\code{coxedExpdur}" and
#' this function will be called by the generic \code{summary} function
#' @param stat Either \code{"mean"} or \code{"median"}
#' @param ... For future methods
#' @details If \code{bootstrap=TRUE} in the call to \code{coxed} then a bootstrapped standard error
#' and confidence interval is reported for the given statistic as well.
#' @return A scalar containing the mean or median duration, or a vector that also includes
#' the bootstrapped standard error and confidence interval for this quantity.
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @export
#' @examples
#' require(survival)
#' mv.surv <- Surv(martinvanberg$formdur, event = rep(1, nrow(martinvanberg)))
#' mv.cox <- coxph(mv.surv ~ postel + prevdef + cont + ident + rgovm + pgovno + tpgovno +
#'      minority, method = "breslow", data = martinvanberg)
#' summary(mv.cox)
#'
#' # NPSF method
#' ed1 <- coxed(mv.cox, method="npsf")
#' ed1$baseline.functions
#' ed1$exp.dur
#' summary(ed1, stat="mean")
#' summary(ed1, stat="median")
#'
#' ed1 <- coxed(mv.cox, method="npsf", bootstrap = TRUE)
#' ed1$exp.dur
#' summary(ed1, stat="mean")
#' summary(ed1, stat="median")
summary.coxedExpdur <- function(object, stat="mean", ...) {
     stopifnot(inherits(object, "coxedExpdur"))
     if(!(stat %in% c("mean", "median"))) stop("stat must be one of mean or median")
     if(stat=="mean") r <- object$mean
     if(stat=="median") r <- object$median
     if(length(r)>1) rownames(r) <- ""
     if(length(r)==1) names(r) <- stat
     r <- round(r, 3)
     return(r)
}
