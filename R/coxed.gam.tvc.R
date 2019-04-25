#' Predict expected durations using the GAM method with time-varying covariates
#'
#' This function is called by \code{\link[coxed]{coxed}} and is not intended to be used by itself.
#' @param cox.model The output from a Cox proportional hazards model estimated
#' with the \code{\link[survival]{coxph}} function in the \code{survival} package
#' or with the \code{\link[rms]{cph}} function in the \code{\link[rms]{rms}}
#' package
#' @param newdata An optional data frame in which to look for variables with
#' which to predict. If omitted, the fitted values are used
#' @param k The number of knots in the GAM smoother. The default is -1, which
#' employs the \code{\link[mgcv]{choose.k}} function from the \code{\link{mgcv}} package
#' to choose the number of knots
#' @param coef A vector of new coefficients to replace the \code{coefficients} attribute
#' of the \code{cox.model}. Used primarily for bootstrapping, to recalculate durations
#' using new coefficients derived from a bootstrapped sample.
#' If \code{NULL}, the original coefficients are employed
#' @param b.ind A vector of observation numbers to pass to the estimation sample to construct
#' the a bootstrapped sample with replacement
#' @param warn If \code{TRUE}, displays warnings, and if \code{FALSE} suppresses them
#' @return Returns a list containing the following components:
#' \tabular{ll}{
#' \code{exp.dur} \tab A vector of predicted mean durations for the estimation sample
#' if \code{newdata} is omitted, or else for the specified new data.  \cr
#' \code{gam.model} \tab Output from the \code{\link[mgcv]{gam}} function in which the durations
#' are fit against the exponentiated linear predictors from the Cox model.\cr
#' \code{gam.data} \tab Fitted values and confidence intervals from the GAM model.\cr
#' }
#' @details This function employs the GAM method of generating expected durations described
#' in Kropko and Harden (2018). See \code{\link[coxed]{coxed.gam}} for details.  This code
#' replicates the code for \code{cox.gam}, but works with data whose structure allows time-varying
#' covariates, and requires using the \code{time2} argument of the \code{\link[survival]{Surv}} function.
#' This function requires the data to be reported as cumulative durations. The GAM model is fit using
#' the ending times for each interval and using only the non-censored observations. Expected durations
#' are drawn from this GAM as with \code{\link[coxed]{coxed.gam}}.
#' @export
#' @seealso \code{\link[mgcv]{gam}}, \code{\link[coxed]{coxed}}, \code{\link[coxed]{coxed.gam}}
#' @references Kropko, J. and Harden, J. J. (2018). Beyond the Hazard Ratio: Generating Expected
#' Durations from the Cox Proportional Hazards Model. \emph{British Journal of Political Science}
#' \url{https://doi.org/10.1017/S000712341700045X}
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @examples
#' bs.surv <- Surv(time = boxsteffensmeier$start, time2 = boxsteffensmeier$te,
#'      event = boxsteffensmeier$cut_hi)
#' bs.cox <- coxph(bs.surv ~ ec + dem + south + iv, data = boxsteffensmeier, method = "breslow")
#' summary(bs.cox)
#'
#' ed <- coxed.gam.tvc(bs.cox)
#' ed$exp.dur
coxed.gam.tvc <- function(cox.model, newdata=NULL, k=-1, coef=NULL, b.ind=NULL, warn=TRUE) {

     exp.xb <- exp(predict(cox.model, type="lp"))

     if(is.null(coef)) gam.data <- data.frame(y=cox.model$y[,2],
                                              failed=cox.model$y[,3],
                                              rank.xb = rank(exp.xb))
     if(!is.null(coef)) gam.data <- data.frame(y=cox.model$y[b.ind,2],
                                               failed=cox.model$y[b.ind,3],
                                               rank.xb = rank(exp.xb[b.ind]))

     gam.model <- gam(y ~ s(rank.xb, bs = "cr", k = k), data = subset(gam.data, failed == 1))

     gam.data <- dplyr::mutate(gam.data,
                        rank.y = rank(y),
                        gam_fit = predict(gam.model, newdata=gam.data),
                        gam_fit_se = predict(gam.model, newdata=gam.data, se.fit=TRUE)$se.fit,
                        gam_fit_95lb = gam_fit + qnorm(.025)*gam_fit_se,
                        gam_fit_95ub = gam_fit + qnorm(.975)*gam_fit_se)

     #New data
     if(!is.null(newdata)){
          exp.xb2 <- exp(predict(cox.model, newdata=newdata, type="lp"))
          rank.xb <- rank.predict(x=exp.xb2, v=exp.xb, warn=warn)
     } else rank.xb <- rank(exp.xb)

     # Generate expected duration from GAM fit
     expect.duration <- predict(gam.model, newdata = data.frame(rank.xb = rank.xb),
                         type = "response")

     return(list(gam.model = gam.model,
                 gam.data = gam.data,
                 exp.dur = expect.duration))

}
