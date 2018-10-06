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
#' @param cluster Cluster variable if bootstrapping is to be done by clusters of
#' observations rather than individual observations. This must be filled in with
#' the case ID if the data are coded with time-varying covariates (using the \code{time2}
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
#' This function rearranges the data so that the time intervals are reported as cumulative durations.
#' For instance, suppose that one observation is observed over three time periods, and thus takes up
#' three rows in the data. On the first row time1=1, time2=4, and event=0; on the second row time1=5,
#' time2=10, and event=0; and on the third row time1=11, time2=13, and event=1. The data are manipulated
#' so that the first row contains duration=3 and event=0, the second row contains duration=10 and event=0,
#' and the third row contains duration=13 and event=1. Then the expected durations are drawn from the Cox model
#' and the NPSF method as with \code{\link[coxed]{coxed.npsf}}.
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
#' ed <- coxed.npsf.tvc(bs.cox, cluster=boxsteffensmeier$caseid)
#' ed$exp.dur
coxed.npsf.tvc <- function(cox.model, newdata=NULL, coef=NULL, b.ind=NULL, cluster) {

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
     if(!is.null(coef)) exp.xb <- exp.xb[b.ind]

     # Compile total failures (only non-censored) at each time point
     if(!is.null(coef)){
          h <- as.data.frame(cbind(y.bs, maxy, failed.bs, exp.xb))
          h <- dplyr::filter(h, y.bs==maxy)
     }
     if(is.null(coef)){
          h <- as.data.frame(cbind(y, maxy, failed, exp.xb))
          h <- dplyr::filter(h, y==maxy)
     }
     h <- dplyr::select(h, -maxy)
     h <- h[order(h[,1]),]
     h <- aggregate(h[,-1], by=list(h[,1]), FUN="sum")
     colnames(h) <- c("time", "total.failures", "exp.xb")

     # Construction of the risk set (includes censored and non-censored observations)
     h[,3] <- rev(cumsum(rev(h[,3])))

     #Construct CBH, baseline survivor and failure CDF
     CBH <- cumsum(h[,2]/h[,3])
     S.bl <- exp(-CBH)
     baseline.functions <- data.frame(time = h$time, cbh = CBH, survivor = S.bl)

     if(!is.null(newdata)){
          exp.xb <- exp(predict(cox.model, newdata=newdata, type="lp"))
     } else{
          exp.xb <- exp(predict(cox.model, type="lp"))
     }

     #Generate EDs for all in-sample observations
     survival <- t(sapply(exp.xb, FUN=function(x){S.bl^x}, simplify=TRUE))
     expect.duration <- apply(survival, 1, FUN=function(x){
          sum(diff(h[,1])*x[-1])
     })
     expect.duration.med <- sum(diff(h[,1])*S.bl^(median(exp.xb))[-1])

     return(list(baseline.functions = baseline.functions,
                 exp.dur = expect.duration))
}
