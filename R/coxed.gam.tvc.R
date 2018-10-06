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
#' @param ties.method A character string specifying how ties are treated,
#' see ‘Details’ in the documentation for \code{\link[base]{rank}}; can be
#' abbreviated
#' @param coef A vector of new coefficients to replace the \code{coefficients} attribute
#' of the \code{cox.model}. Used primarily for bootstrapping, to recalculate durations
#' using new coefficients derived from a bootstrapped sample.
#' If \code{NULL}, the original coefficients are employed
#' @param b.ind A vector of observation numbers to pass to the estimation sample to construct
#' the a bootstrapped sample with replacement
#' @param warn If \code{TRUE}, displays warnings, and if \code{FALSE} suppresses them
#' @param cluster Cluster variable if bootstrapping is to be done by clusters of
#' observations rather than individual observations. This must be filled in with
#' the case ID if the data are coded with time-varying covariates (using the \code{time2}
#' argument in the \code{\link[survival]{Surv}} function)
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
#' This function rearranges the data so that the time intervals are reported as cumulative durations.
#' For instance, suppose that one observation is observed over three time periods, and thus takes up
#' three rows in the data. On the first row time1=1, time2=4, and event=0; on the second row time1=5,
#' time2=10, and event=0; and on the third row time1=11, time2=13, and event=1. The data are manipulated
#' so that the first row contains duration=3 and event=0, the second row contains duration=10 and event=0,
#' and the third row contains duration=13 and event=1. Then the GAM model is fit using the non-censored
#' observations, and expected durations are drawn from this GAM as with \code{\link[coxed]{coxed.gam}}.
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
#' ed <- coxed.gam.tvc(bs.cox, cluster=boxsteffensmeier$caseid)
#' ed$exp.dur
coxed.gam.tvc <- function(cox.model, newdata=NULL, k=-1, ties.method="random",
                              coef=NULL, b.ind=NULL, warn=TRUE, cluster) {

     if(!is.null(coef)){
          df <- data.frame(y = cox.model$y[b.ind,2], id = cluster[b.ind])
          df <- dplyr::group_by(df, id)
          df <- dplyr::mutate(df, maxy = max(y))
          df <- dplyr::ungroup(df)
          y.bs <- df$y
          maxy <- df$maxy
          failed.bs <- cox.model$y[b.ind,3]
          cox.model$coefficients <- coef
     }
     df <- data.frame(y = cox.model$y[,2], id = cluster)
     df <- dplyr::group_by(df, id)
     df <- dplyr::mutate(df, maxy = max(y))
     df <- dplyr::ungroup(df)
     y <- df$y
     maxy <- df$maxy
     failed <- cox.model$y[,3]
     exp.xb <- exp(predict(cox.model, type="lp"))
     if(!is.null(coef)) exp.xb.bs <- exp.xb[b.ind]

     gam.data <- data.frame(y=y, failed=failed, rank.xb = rank(exp.xb, ties.method = ties.method))

     if(!is.null(coef)){
          gam.data.bs <- data.frame(y=y.bs, failed=failed.bs, rank.xb = rank(exp.xb.bs, ties.method = ties.method))
          gam.model <- gam(y ~ s(rank.xb, bs = "cr", k = k), data = subset(gam.data.bs, failed == 1 & y.bs==maxy))
     } else {
          gam.model <- gam(y ~ s(rank.xb, bs = "cr", k = k), data = subset(gam.data, failed == 1 & y==maxy))
     }

     gam.data <- dplyr::mutate(gam.data,
                        rank.y = rank(y, ties.method = ties.method),
                        gam_fit = predict(gam.model, newdata=gam.data),
                        gam_fit_se = predict(gam.model, newdata=gam.data, se.fit=TRUE)$se.fit,
                        gam_fit_95lb = gam_fit + qnorm(.025)*gam_fit_se,
                        gam_fit_95ub = gam_fit + qnorm(.975)*gam_fit_se)

     #New data
     if(!is.null(newdata)){
          exp.xb2 <- exp(predict(cox.model, newdata=newdata, type="lp"))
          rank.xb <- rank.predict(x=exp.xb2, v=exp.xb, ties.method=ties.method, warn=warn)
     } else rank.xb <- rank(exp.xb, ties.method = ties.method)

     # Generate expected duration from GAM fit
     expect.duration <- predict(gam.model, newdata = data.frame(rank.xb = rank.xb),
                         type = "response")

     return(list(gam.model = gam.model,
                 gam.data = gam.data,
                 exp.dur = expect.duration))

}
