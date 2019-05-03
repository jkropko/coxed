#' Predict expected durations using the GAM method with time-varying covariates
#'
#' This function is called by \code{\link[coxed]{coxed}} and is not intended to be used by itself.
#' @param cox.model The output from a Cox proportional hazards model estimated
#' with the \code{\link[survival]{coxph}} function in the \code{survival} package
#' or with the \code{\link[rms]{cph}} function in the \code{\link[rms]{rms}}
#' package
#' @param newdata An optional data frame in which to look for variables with
#' which to predict. If omitted, the fitted values are used
#' @param coef A vector of new coefficients to replace the \code{coefficients} attribute
#' of the \code{cox.model}. Used primarily for bootstrapping, to recalculate durations
#' using new coefficients derived from a bootstrapped sample.
#' If \code{NULL}, the original coefficients are employed
#' @param b.ind A vector of observation numbers to pass to the estimation sample to construct
#' the a bootstrapped sample with replacement
#' @return Returns a list containing the following components:
#' \tabular{ll}{
#' \code{exp.dur} \tab A vector of predicted mean durations for the estimation sample
#' if \code{newdata} is omitted, or else for the specified new data. \cr
#' \code{baseline.functions} \tab The estimated cumulative baseline hazard function and survivor function. \cr
#' }
#' @details This function employs the NPSF method of generating expected durations described
#' in Kropko and Harden (2018). See \code{\link[coxed]{coxed.npsf}} for details.  This code
#' replicates the code for \code{cox.npsf}, but works with data whose structure allows time-varying
#' covariates, and requires using the \code{time2} argument of the \code{\link[survival]{Surv}} function.
#' This function requires the data to be reported as cumulative durations. The cumulative baseline hazard
#' function model is estimated using the ending times for each interval. Then the expected durations are
#' drawn from the Cox model and the NPSF method as with \code{\link[coxed]{coxed.npsf}}.
#' @seealso \code{\link[coxed]{coxed}}, \code{\link[coxed]{coxed.npsf}}
#' @references Kropko, J. and Harden, J. J. (2018). Beyond the Hazard Ratio: Generating Expected Durations
#' from the Cox Proportional Hazards Model. \emph{British Journal of Political Science}
#' \url{https://doi.org/10.1017/S000712341700045X}
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @export
#' @examples
#' bs.surv <- Surv(time = boxsteffensmeier$start, time2 = boxsteffensmeier$te,
#'      event = boxsteffensmeier$cut_hi)
#' bs.cox <- coxph(bs.surv ~ ec + dem + south + iv, data = boxsteffensmeier, method = "breslow")
#'
#' ed <- coxed.npsf.tvc(bs.cox)
#' ed$exp.dur
coxed.npsf.tvc <- function(cox.model, newdata=NULL, coef=NULL, b.ind=NULL) {

     start <- ceiling(cox.model$y[,1]) ## ceiling() to deal with decimals
     end <- ceiling(cox.model$y[,2])
     failed <- cox.model$y[,3]
     exp.xb <- exp(predict(cox.model, type="lp"))
     if(!is.null(coef)){
          start <- start[b.ind]
          end <- end[b.ind]
          failed <- failed[b.ind]
          exp.xb <- exp.xb[b.ind]
     }

     h <- as.data.frame(cbind(start, end, failed, exp.xb))
     diff <- h$end - h$start
     h <- h[rep(1:nrow(h), diff),]
     h$time <- h$start + sequence(diff)
     h$failed <- ifelse(h$time==h$end, h$failed, 0)

     h <- dplyr::group_by(h, time)
     h <- dplyr::summarize(h, d = sum(failed),
                           exp.xb = sum(exp.xb))
     #Construct CBH, baseline survivor and failure CDF
     CBH <- cumsum(h$d / h$exp.xb)
     S.bl <- exp(-CBH)
     baseline.functions <- data.frame(time = h$time, cbh = CBH, survivor = S.bl)

     if(!is.null(newdata)) exp.xb <- exp(predict(cox.model, newdata=newdata, type="lp"))

     #Generate EDs for all observations in exp.xb
     expect.duration <- sapply(exp.xb, FUN=function(x){
          sum(S.bl^x)
     })

     return(list(baseline.functions = baseline.functions,
                 exp.dur = expect.duration))
}
