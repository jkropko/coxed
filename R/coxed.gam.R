#' Predict expected durations using the GAM method
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
#' if \code{newdata} is omitted, or else for the specified new data. \cr
#' \code{gam.model} \tab Output from the \code{\link[mgcv]{gam}} function in which the durations
#' are fit against the exponentiated linear predictors from the Cox model.\cr
#' \code{gam.data} \tab Fitted values and confidence intervals from the GAM model.\cr
#' }
#' @details This function employs the GAM method of generating expected durations described
#' in Kropko and Harden (2018), which proceeds according to five steps. First, it uses coefficient
#' estimates from the Cox model, so researchers must first estimate the model just as they always have.
#' Then the method computes expected values of risk for each observation by matrix-multiplying the
#' covariates by the estimated coefficients from the model, then exponentiating the result. This creates
#' the exponentiated linear predictor (ELP). Then the observations are ranked from smallest to largest
#' according to their values of the ELP. This ranking is interpreted as the expected order of failure;
#' the larger the value of the ELP, the sooner the model expects that observation to fail, relative to
#' the other observations.
#'
#' The next step is to connect the model's expected risk for each observation (ELP) to duration time
#' (the observed durations). A \code{\link[mgcv]{gam}} fits a model to data by using a series of locally-estimated polynomial
#' splines set by the user (see, for example, Wood, Pya, and Saefken 2016). It is a flexible means of allowing for
#' the possibility of nonlinear relationships between variables. \code{coxed.gam} uses a GAM to model the observed
#' utilizes a cubic regression spline to draw a smoothed line summarizing the bivariate relationship between
#' the observed durations and the ranks. The GAM fit can be used directly to compute expected durations, given the covariates, for each observation
#' in the data.
#' @export
#' @seealso \code{\link[mgcv]{gam}}, \code{\link[coxed]{coxed}}
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @references Kropko, J. and Harden, J. J. (2018). Beyond the Hazard Ratio: Generating Expected
#' Durations from the Cox Proportional Hazards Model. \emph{British Journal of Political Science}
#' \url{https://doi.org/10.1017/S000712341700045X}
#'
#' Wood, S.N., N. Pya and B. Saefken (2016). Smoothing parameter and model selection for general smooth models (with discussion). \emph{Journal of the American Statistical Association} \strong{111}, 1548-1575
#' \url{http://dx.doi.org/10.1080/01621459.2016.1180986}
#' @examples
#' mv.surv <- Surv(martinvanberg$formdur, event = rep(1, nrow(martinvanberg)))
#' mv.cox <- coxph(mv.surv ~ postel + prevdef + cont + ident + rgovm +
#' pgovno + tpgovno + minority, method = "breslow", data = martinvanberg)
#'
#' ed <- coxed.gam(mv.cox)
#' summary(ed$gam.data)
#' summary(ed$gam.model)
#' ed$exp.dur
#'
#' #Plotting the GAM fit
#' \dontrun{require(ggplot2)
#' ggplot(ed$gam.data, aes(x=rank.xb, y=y)) +
#'     geom_point() +
#'     geom_line(aes(x=rank.xb, y=gam_fit)) +
#'     geom_ribbon(aes(ymin=gam_fit_95lb, ymax=gam_fit_95ub), alpha=.5) +
#'     xlab("Cox model LP rank (smallest to largest)") +
#'     ylab("Duration")
#' }
#'
#' #Running coxed.gam() on a bootstrap sample and with new coefficients
#' bsample <- sample(1:nrow(martinvanberg), nrow(martinvanberg), replace=TRUE)
#' newcoefs <- rnorm(8)
#' ed2 <- coxed.gam(mv.cox, b.ind=bsample, coef=newcoefs)
coxed.gam <- function(cox.model, newdata=NULL, k=-1, coef=NULL, b.ind=NULL, warn=TRUE) {

     if(!is.null(coef)){
          y.bs <- cox.model$y[b.ind,1]
          failed.bs <- cox.model$y[b.ind,2]
          cox.model$coefficients <- coef
     }
     y <- cox.model$y[,1]
     failed <- cox.model$y[,2]
     exp.xb <- exp(predict(cox.model, type="lp"))
     if(!is.null(coef)) exp.xb.bs <- exp.xb[b.ind]

     gam.data <- data.frame(y=y, failed=failed, rank.xb = rank(exp.xb))

     if(!is.null(coef)){
          gam.data.bs <- data.frame(y=y.bs, failed=failed.bs, rank.xb = rank(exp.xb.bs))
          gam.model <- gam(y ~ s(rank.xb, bs = "cr", k = k), data = subset(gam.data.bs, failed == 1))
     } else {
          gam.model <- gam(y ~ s(rank.xb, bs = "cr", k = k), data = subset(gam.data, failed == 1))
     }

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
