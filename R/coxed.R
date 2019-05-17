#' Expected durations and marginal changes in expected duration from the
#' Cox proportional hazards model
#'
#' \code{coxed()} returns expected durations for every observation in the data
#' used to fit the model, or in new data, or returns the mean or median of these
#' durations, or differences in duration for two pre-defined covariate profiles.
#' Standard errors and confidence intervals for all quantities produced by
#' \code{coxed()} are calculated via bootstrapping.
#'
#' @param cox.model The output from a Cox proportional hazards model estimated
#' with the \code{\link[survival]{coxph}} function in the \code{survival} package
#' or with the \code{\link[rms]{cph}} function in the \code{\link[rms]{rms}}
#' package
#' @param newdata An optional data frame in which to look for variables with
#' which to predict. If omitted, the fitted values are used
#' @param newdata2 An optional data frame that can only be specified if
#' \code{newdata} is not omitted, and must have the same dimensions as \code{newdata}.
#' If specified, marginal changes are calculated by subtracting the expected
#' durations for \code{newdata2} from the expected durations for \code{newdata}
#' @param bootstrap Should bootstrapped standard errors and confidence intervals
#' be calculated?
#' @param method If "npsf" (the default), expected durations are calculated using
#' the non-parametric step function approach described in Kropko and Harden
#' (2018). If "gam", expected durations are calculated using the GAM method
#' @param k The number of knots in the GAM smoother. The default is -1, which
#' employs the \code{\link[mgcv]{choose.k}} function from the \code{\link{mgcv}} package
#' to choose the number of knots
#' @param B Number of bootstrap simulation iterations
#' @param confidence  If "studentized" (the default), bootstrapped CIs
#' are calculated from the tails of a normal distribution where the mean and
#' standard deviation are the point estimate and boostrapped SE of each duration
#' estimate. If "empirical", bootstrapped confidence intervals are calculated empirically.
#' If "bca", bootstrapped confidence intervals are calculated using the bias-correction
#' and acceleration method described by DiCiccio and Efron (1996).
#' @param level The level of the confidence interval to calculate (default is
#' .95 for a 95 percent confidence interval)
#' @param id Cluster variable if bootstrapping is to be done by clusters of
#' observations rather than individual observations. If the data are coded
#' with time-varying covariates (using the \code{time2} argument in the
#' \code{\link[survival]{Surv}} function), this variable must be the ID variable
#' in the \emph{data that are used to estimate the Cox PH model}, and not the ID
#' variable in new data.
#' @param ... Additional arguments to be passed to the \code{\link[coxed]{bootcov2}}
#' function, an adaptation of the \code{\link[rms]{bootcov}} function in the
#' \code{\link{rms}} package
#'
#' @details The \code{coxed} function generates expected durations for individual
#' observations and/or marginal changes in expected duration given a change in a covariate
#' from the Cox proportional hazards model. Specifically, the methods can compute (1) the
#' expected duration for each observation used to fit the Cox model, given the covariates,
#' (2) the expected duration for a "new" observation with a covariate profile set by the
#' analyst, or (3) the first difference, or change, in expected duration given two new data frames.
#'
#' There are two different methods, described in Kropko and Harden (2018), of generating duration-based quantities in the package.
#' The first method calculates expected durations by using a nonparametric estimate of the
#' baseline hazard and survivor functions (see \code{\link[coxed]{coxed.npsf}} for details).
#' The second method employs a generalized additive model
#' (GAM) to map the model's estimated linear predictor values to duration times (see \code{\link[coxed]{coxed.gam}} for details). Both
#' methods are also implemented for data structures with time-varying covariates
#' (see \code{\link[coxed]{coxed.npsf.tvc}} and \code{\link[coxed]{coxed.gam.tvc}}).
#' @return \code{coxed} returns an object of \code{\link[base]{class}} "coxedExpdur" or
#' "coxedMargin", which is a list containing some of the following components, depending
#' on the implementation of \code{coxed}:
#' \tabular{ll}{
#' \code{exp.dur} \tab A vector of predicted mean durations for the estimation sample
#' if \code{newdata} is omitted, or else for the specified new data.  If \code{bootstrap}
#' is \code{TRUE} bootstrapped standard errors are also provided, as well as the confidence
#' interval requested by \code{level}. \cr
#' \code{mean} \tab The mean of the predicted durations. If \code{bootstrap}
#' is \code{TRUE} bootstrapped standard errors are also provided, as well as the confidence
#' interval requested by \code{level}. \cr
#' \code{median} \tab The median of the predicted durations. If \code{bootstrap}
#' is \code{TRUE} bootstrapped standard errors are also provided, as well as the confidence
#' interval requested by \code{level}. \cr
#' \code{baseline.functions} \tab The estimated cumulative baseline hazard function and survivor function. \cr
#' \code{gam.model} \tab Output from the \code{\link[mgcv]{gam}} function in which the durations
#' are fit against the exponentiated linear predictors from the Cox model.\cr
#' \code{gam.data} \tab Fitted values and confidence intervals from the GAM model.\cr
#' \code{exp.dur1} \tab  A vector of predicted mean durations for the observations in \code{newdata1}
#' when calculating marginal effects. \cr
#' \code{exp.dur2} \tab A vector of predicted mean durations for the observations in \code{newdata2}
#' when calculating marginal effects. \cr
#' \code{mean1} \tab The mean of the predicted mean durations for the observations in \code{newdata1}
#' when calculating marginal effects. \cr
#' \code{mean2} \tab The mean of the predicted mean durations for the observations in \code{newdata2}
#' when calculating marginal effects. \cr
#' \code{median1} \tab The median of the predicted mean durations for the observations in \code{newdata1}
#' when calculating marginal effects. \cr
#' \code{median2} \tab The median of the predicted mean durations for the observations in \code{newdata2}
#' when calculating marginal effects. \cr
#' \code{diff} \tab A vector of the difference between the predicted mean durations for each
#' observation under the covariate profile in \code{newdata2} and the covariate profile in \code{newdata1}.\cr
#' \code{mean.diff} \tab The mean of the differences in duration across observations. \cr
#' \code{median.diff} \tab The median of the differences in duration across observations. \cr
#' }
#' @references Kropko, J. and Harden, J. J. (2018). Beyond the Hazard Ratio: Generating Expected Durations from the
#' Cox Proportional Hazards Model. \emph{British Journal of Political Science} \url{https://doi.org/10.1017/S000712341700045X}
#'
#' DiCiccio, T. J. and B. Efron. (1996). Bootstrap Confidence Intervals. \emph{Statistical Science}.
#' 11(3): 189–212. \url{https://doi.org/10.1214/ss/1032280214}
#' @seealso \code{\link[survival]{coxph}}, \code{\link[rms]{cph}}, \code{\link[coxed]{bootcov2}},
#' \code{\link[coxed]{coxed.gam}}, \code{\link[coxed]{coxed.gam.tvc}}, \code{\link[coxed]{coxed.npsf}},
#' \code{\link[coxed]{coxed.npsf.tvc}}
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @export
#' @examples
#' mv.surv <- Surv(martinvanberg$formdur, event = rep(1, nrow(martinvanberg)))
#' mv.cox <- coxph(mv.surv ~ postel + prevdef + cont + ident + rgovm + pgovno +
#'      tpgovno + minority, method = "breslow", data = martinvanberg)
#' summary(mv.cox)
#'
#' # NPSF method
#' ed1 <- coxed(mv.cox, method="npsf")
#' ed1$baseline.functions
#' ed1$exp.dur
#' summary(ed1, stat="mean")
#' summary(ed1, stat="median")
#'
#' \dontrun{ed1 <- coxed(mv.cox, method="npsf", bootstrap = TRUE)
#' ed1$exp.dur
#' summary(ed1, stat="mean")
#' summary(ed1, stat="median")
#' }
#'
#' me <- coxed(mv.cox, method="npsf", bootstrap = FALSE,
#'             newdata = dplyr::mutate(martinvanberg, pgovno=1),
#'             newdata2 = dplyr::mutate(martinvanberg, pgovno=6))
#' summary(me, stat="mean")
#'
#' # GAM method
#' ed2 <- coxed(mv.cox, method="gam")
#' summary(ed2$gam.data)
#' summary(ed2$gam.model)
#' ed2$exp.dur
#' summary(ed2, stat="mean")
#'
#' \dontrun{me <- coxed(mv.cox, method="gam", bootstrap = TRUE,
#'             newdata = dplyr::mutate(martinvanberg, pgovno=1),
#'             newdata2 = dplyr::mutate(martinvanberg, pgovno=6))
#' summary(me, stat="mean")
#' summary(me, stat="median")
#' }
#'
#' #Plotting the GAM fit
#' \dontrun{ggplot(ed2$gam.data, aes(x=rank.xb, y=y)) +
#'     geom_point() +
#'     geom_line(aes(x=rank.xb, y=gam_fit)) +
#'     geom_ribbon(aes(ymin=gam_fit_95lb, ymax=gam_fit_95ub), alpha=.5) +
#'     xlab("Cox model LP rank (smallest to largest)") +
#'     ylab("Duration")
#' }
#'
#' #Time-varying covariates
#' bs.surv <- Surv(time = boxsteffensmeier$start, time2 = boxsteffensmeier$te,
#'      event = boxsteffensmeier$cut_hi)
#' bs.cox <- coxph(bs.surv ~ ec + dem + south + iv, data = boxsteffensmeier, method = "breslow")
#' summary(bs.cox)
#'
#' ed1 <- coxed(bs.cox, method="npsf", id=boxsteffensmeier$caseid)
#' ed1$exp.dur
#' summary(ed1, stat="mean")

coxed <- function(cox.model, newdata=NULL, newdata2=NULL, bootstrap=FALSE, method="npsf",
                  k=-1, B = 200, confidence="studentized",
                  level=.95, id=NULL, ...){

     tvc <- (ncol(cox.model$y)==3)
     if(is.null(newdata) & !is.null(newdata2)) stop("newdata2 must only be specified if newdata is also specified")
     if(!method %in% c("gam", "npsf")) stop("method must be one of 'gam' and 'npsf'")
     if(!is.null(newdata2) & !all(dim(newdata)==dim(newdata2))) stop("newdata and newdata2 must have the same dimensions")
     if(!confidence %in% c("studentized", "empirical", "bca")) stop("confidence must be one of 'studentized', 'empirical', or 'bca'")
     if(tvc & is.null(id)) stop("id must be filled in with the case ID in the data used to estimate the Cox PH model if you have time-varying covariates.")

     #First data frame (newdata), or estimation sample if NULL
     if(method=="gam"){
          if(!tvc) dur1 <- coxed.gam(cox.model, newdata=newdata, k=k)
          if(tvc) dur1 <- coxed.gam.tvc(cox.model, newdata=newdata, k=k)
          dur1.pe <- dur1$exp.dur
          gam.model <- dur1$gam.model
          gam.data <- dur1$gam.data
          m1 <- mean(gam.data$y, na.rm=TRUE)
          m2 <- mean(gam.data$y[gam.data$failed==1], na.rm=TRUE)
          if(m1/m2 > 1.25) warning("Inclusion of censored cases increases mean duration by more than 25%. Consider using method='npsf' instead")
     }
     if(method=="npsf"){
          if(!tvc) dur1 <- coxed.npsf(cox.model, newdata=newdata)
          if(tvc) dur1 <- coxed.npsf.tvc(cox.model, newdata=newdata)
          dur1.pe <- dur1$exp.dur
          baseline.functions <- dur1$baseline.functions
     }
     exp.dur <- data.frame(exp.dur=dur1.pe)

     #Second data frame (newdata2)
     if(!is.null(newdata2)){
          if(method=="gam" & !tvc) dur2.pe <- coxed.gam(cox.model, newdata=newdata2,
                                                        k=k)$exp.dur
          if(method=="gam" & tvc) dur2.pe <- coxed.gam.tvc(cox.model, newdata=newdata2,
                                                           k=k)$exp.dur
          if(method=="npsf" & !tvc) dur2.pe <- coxed.npsf(cox.model, newdata=newdata2)$exp.dur
          if(method=="npsf" & tvc) dur2.pe <- coxed.npsf.tvc(cox.model, newdata=newdata2)$exp.dur
          diff <- dur2.pe - dur1.pe
          exp.dur <- data.frame(exp.dur1 = dur1.pe, exp.dur2=dur2.pe, difference=diff)
     }

     if(!bootstrap){
          res <- list(exp.dur = exp.dur, mean = colMeans(exp.dur, na.rm=TRUE),
                      median = apply(exp.dur, 2, FUN=function(x){
                           median(x, na.rm=TRUE)
                      }))
     }

     if(bootstrap){
          if(!tvc) boot.model <- bootcov2(rms::cph(Surv(as.numeric(cox.model$y[ , 1]),
                                                        as.numeric(cox.model$y[ , 2]),
                                                        type = "right") ~
                                                        model.matrix(cox.model),
                                                   x = TRUE, y = TRUE), B = B, ...)
          if(tvc){
               id <- id[as.numeric(rownames(model.frame(cox.model)))]
               boot.cph <- rms::cph(Surv(as.numeric(cox.model$y[ , 1]),
                                         as.numeric(cox.model$y[ , 2]),
                                         as.numeric(cox.model$y[ , 3]),
                                         type = "counting") ~
                                         model.matrix(cox.model),
                                    x = TRUE, y = TRUE)
               class(boot.cph) <- "tvc"
               boot.model <- bootcov2(boot.cph, B = B, cluster=id, ...)
          }
          bs.coef <- boot.model$boot.Coef
          bs.obs <- na.omit(boot.model$b.ind)

          if(is.null(newdata2)){
               exp.dur.mat <- matrix(numeric(), nrow(exp.dur), nrow(bs.coef))
               for(i in 1:nrow(bs.coef)){
                    if(method=="gam" & !tvc) dur1.pe <- coxed.gam(cox.model, newdata=newdata, warn=FALSE,
                                                                  k=k,
                                                                  coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                    if(method=="gam" & tvc) dur1.pe <- coxed.gam.tvc(cox.model, newdata=newdata, warn=FALSE,
                                                                     k=k,
                                                                     coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                    if(method=="npsf" & !tvc) dur1.pe <- coxed.npsf(cox.model, newdata=newdata,
                                                                    coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                    if(method=="npsf" & tvc) dur1.pe <- coxed.npsf.tvc(cox.model, newdata=newdata,
                                                                       coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                    exp.dur.mat[,i] <- dur1.pe
               }
               exp.dur.se <- apply(exp.dur.mat, 1, FUN=function(x){
                    sd(x, na.rm=TRUE)
               })
               mean.vec <- colMeans(exp.dur.mat, na.rm=TRUE)
               mean.se <- sd(mean.vec, na.rm=TRUE)
               median.vec <- apply(exp.dur.mat, 2, FUN=function(x){
                    median(x, na.rm=TRUE)
               })
               median.se <- sd(median.vec, na.rm=TRUE)
               if(confidence=="empirical"){
                    exp.dur.lb <- apply(exp.dur.mat, 1, FUN=function(x){
                         quantile(x, (1-level)/2, names=FALSE, na.rm=TRUE)
                    })
                    exp.dur.ub <- apply(exp.dur.mat, 1, FUN=function(x){
                         quantile(x, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    })
                    mean.lb <- quantile(mean.vec, (1-level)/2, names=FALSE, na.rm=TRUE)
                    mean.ub <- quantile(mean.vec, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    median.lb <- quantile(median.vec, (1-level)/2, names=FALSE, na.rm=TRUE)
                    median.ub <- quantile(median.vec, level + (1-level)/2, names=FALSE, na.rm=TRUE)
               } else if(confidence=="studentized"){
                    exp.dur.lb <- exp.dur$exp.dur + qnorm((1-level)/2)*exp.dur.se
                    exp.dur.ub <- exp.dur$exp.dur + qnorm(level + (1-level)/2)*exp.dur.se
                    mean.lb <- mean(exp.dur$exp.dur, na.rm=TRUE) + qnorm((1-level)/2)*mean.se
                    mean.ub <- mean(exp.dur$exp.dur, na.rm=TRUE) + qnorm(level + (1-level)/2)*mean.se
                    median.lb <- median(exp.dur$exp.dur, na.rm=TRUE) + qnorm((1-level)/2)*median.se
                    median.ub <- median(exp.dur$exp.dur, na.rm=TRUE) + qnorm(level + (1-level)/2)*median.se
               } else if(confidence=="bca"){
                    exp.dur.bca <- apply(exp.dur.mat, 1, FUN=function(x){
                         bca(x, conf.level = level)
                    })
                    exp.dur.lb <- exp.dur.bca[1,]
                    exp.dur.ub <- exp.dur.bca[2,]
                    mean.bca <- bca(mean.vec, conf.level = level)
                    mean.lb <- mean.bca[1]
                    mean.ub <- mean.bca[2]
                    median.bca <- bca(median.vec, conf.level = level)
                    median.lb <- median.bca[1]
                    median.ub <- median.bca[2]
               }
               res <- list(exp.dur = data.frame(exp.dur = exp.dur$exp.dur, bootstrap.se = exp.dur.se,
                                                lb = exp.dur.lb, ub = exp.dur.ub),
                           mean = data.frame(mean = mean(exp.dur$exp.dur, na.rm=TRUE), bootstrap.se = mean.se,
                                             lb = mean.lb, ub = mean.ub),
                           median = data.frame(median = median(exp.dur$exp.dur, na.rm=TRUE), bootstrap.se = median.se,
                                               lb = median.lb, ub = median.ub))
          } else {
               exp.dur1.mat <- matrix(numeric(), nrow(exp.dur), nrow(bs.coef))
               exp.dur2.mat <- matrix(numeric(), nrow(exp.dur), nrow(bs.coef))
               diff.mat <- matrix(numeric(), nrow(exp.dur), nrow(bs.coef))
               for(i in 1:nrow(bs.coef)){
                    if(method=="gam" & !tvc){
                         warn <- (i==1)
                         dur1.pe <- coxed.gam(cox.model, newdata=newdata,
                                              k=k, warn=warn,
                                              coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                         dur2.pe <- coxed.gam(cox.model, newdata=newdata2,
                                              k=k, warn=warn,
                                              coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                    }
                    if(method=="gam" & tvc){
                         warn <- (i==1)
                         dur1.pe <- coxed.gam.tvc(cox.model, newdata=newdata,
                                                  k=k, warn=warn,
                                                  coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                         dur2.pe <- coxed.gam.tvc(cox.model, newdata=newdata2,
                                                  k=k, warn=warn,
                                                  coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                    }
                    if(method=="npsf" & !tvc){
                         dur1.pe <- coxed.npsf(cox.model, newdata=newdata,
                                               coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                         dur2.pe <- coxed.npsf(cox.model, newdata=newdata2,
                                               coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                    }
                    if(method=="npsf" & tvc){
                         dur1.pe <- coxed.npsf.tvc(cox.model, newdata=newdata,
                                                   coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                         dur2.pe <- coxed.npsf.tvc(cox.model, newdata=newdata2,
                                                   coef=bs.coef[i,], b.ind=bs.obs[,i])$exp.dur
                    }
                    exp.dur1.mat[,i] <- dur1.pe
                    exp.dur2.mat[,i] <- dur2.pe
                    diff.mat[,i] <- dur2.pe - dur1.pe
               }
               exp.dur1.se <- apply(exp.dur1.mat, 1, FUN=function(x){
                    sd(x, na.rm=TRUE)
               })
               mean1.vec <- colMeans(exp.dur1.mat, na.rm=TRUE)
               mean1.se <- sd(mean1.vec, na.rm=TRUE)
               median1.vec <- apply(exp.dur1.mat, 2, FUN=function(x){
                    median(x, na.rm=TRUE)
               })
               median1.se <- sd(median1.vec, na.rm=TRUE)
               exp.dur2.se <- apply(exp.dur2.mat, 1, FUN=function(x){
                    sd(x, na.rm=TRUE)
               })
               mean2.vec <- colMeans(exp.dur2.mat, na.rm=TRUE)
               mean2.se <- sd(mean2.vec, na.rm=TRUE)
               median2.vec <- apply(exp.dur2.mat, 2, FUN=function(x){
                    median(x, na.rm=TRUE)
               })
               median2.se <- sd(median2.vec, na.rm=TRUE)
               diff.se <- apply(diff.mat, 1, FUN=function(x){
                    sd(x, na.rm=TRUE)
               })
               mean.diff.vec <- colMeans(diff.mat, na.rm=TRUE)
               mean.diff.se <- sd(mean.diff.vec, na.rm=TRUE)
               median.diff.vec <- apply(diff.mat, 2, FUN=function(x){
                    median(x, na.rm=TRUE)
               })
               median.diff.se <- sd(median.diff.vec, na.rm=TRUE)
               if(confidence=="empirical"){
                    exp.dur1.lb <- apply(exp.dur1.mat, 1, FUN=function(x){
                         quantile(x, (1-level)/2, names=FALSE, na.rm=TRUE)
                    })
                    exp.dur1.ub <- apply(exp.dur1.mat, 1, FUN=function(x){
                         quantile(x, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    })
                    mean1.lb <- quantile(mean1.vec, (1-level)/2, names=FALSE, na.rm=TRUE)
                    mean1.ub <- quantile(mean1.vec, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    median1.lb <- quantile(median1.vec, (1-level)/2, names=FALSE, na.rm=TRUE)
                    median1.ub <- quantile(median1.vec, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    exp.dur2.lb <- apply(exp.dur2.mat, 1, FUN=function(x){
                         quantile(x, (1-level)/2, names=FALSE, na.rm=TRUE)
                    })
                    exp.dur2.ub <- apply(exp.dur2.mat, 1, FUN=function(x){
                         quantile(x, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    })
                    mean2.lb <- quantile(mean2.vec, (1-level)/2, names=FALSE, na.rm=TRUE)
                    mean2.ub <- quantile(mean2.vec, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    median2.lb <- quantile(median2.vec, (1-level)/2, names=FALSE, na.rm=TRUE)
                    median2.ub <- quantile(median2.vec, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    diff.lb <- apply(diff.mat, 1, FUN=function(x){
                         quantile(x, (1-level)/2, names=FALSE, na.rm=TRUE)
                    })
                    diff.ub <- apply(diff.mat, 1, FUN=function(x){
                         quantile(x, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    })
                    mean.diff.lb <- quantile(mean.diff.vec, (1-level)/2, names=FALSE, na.rm=TRUE)
                    mean.diff.ub <- quantile(mean.diff.vec, level + (1-level)/2, names=FALSE, na.rm=TRUE)
                    median.diff.lb <- quantile(median.diff.vec, (1-level)/2, names=FALSE, na.rm=TRUE)
                    median.diff.ub <- quantile(median.diff.vec, level + (1-level)/2, names=FALSE, na.rm=TRUE)
               } else if(confidence=="studentized"){
                    exp.dur1.lb <- exp.dur[,1] + qnorm((1-level)/2)*exp.dur1.se
                    exp.dur1.ub <- exp.dur[,1] + qnorm(level + (1-level)/2)*exp.dur1.se
                    mean1.lb <- mean(exp.dur[,1], na.rm=TRUE) + qnorm((1-level)/2)*mean1.se
                    mean1.ub <- mean(exp.dur[,1], na.rm=TRUE) + qnorm(level + (1-level)/2)*mean1.se
                    median1.lb <- median(exp.dur[,1], na.rm=TRUE) + qnorm((1-level)/2)*median1.se
                    median1.ub <- median(exp.dur[,1], na.rm=TRUE) + qnorm(level + (1-level)/2)*median1.se
                    exp.dur2.lb <- exp.dur[,2] + qnorm((1-level)/2)*exp.dur2.se
                    exp.dur2.ub <- exp.dur[,2] + qnorm(level + (1-level)/2)*exp.dur2.se
                    mean2.lb <- mean(exp.dur[,2], na.rm=TRUE) + qnorm((1-level)/2)*mean2.se
                    mean2.ub <- mean(exp.dur[,2], na.rm=TRUE) + qnorm(level + (1-level)/2)*mean2.se
                    median2.lb <- median(exp.dur[,2], na.rm=TRUE) + qnorm((1-level)/2)*median2.se
                    median2.ub <- median(exp.dur[,2], na.rm=TRUE) + qnorm(level + (1-level)/2)*median2.se
                    diff.lb <- exp.dur[,3] + qnorm((1-level)/2)*diff.se
                    diff.ub <- exp.dur[,3] + qnorm(level + (1-level)/2)*diff.se
                    mean.diff.lb <- mean(exp.dur[,3], na.rm=TRUE) + qnorm((1-level)/2)*mean.diff.se
                    mean.diff.ub <- mean(exp.dur[,3], na.rm=TRUE) + qnorm(level + (1-level)/2)*mean.diff.se
                    median.diff.lb <- median(exp.dur[,3], na.rm=TRUE) + qnorm((1-level)/2)*median.diff.se
                    median.diff.ub <- median(exp.dur[,3], na.rm=TRUE) + qnorm(level + (1-level)/2)*median.diff.se
               } else if(confidence=="bca"){
                    exp.dur1.bca <- apply(exp.dur1.mat, 1, FUN=function(x){
                         bca(x, conf.level = level)
                    })
                    exp.dur1.lb <- exp.dur1.bca[1,]
                    exp.dur1.ub <- exp.dur1.bca[2,]
                    mean1.bca <- bca(mean1.vec, conf.level = level)
                    mean1.lb <- mean1.bca[1]
                    mean1.ub <- mean1.bca[2]
                    median1.bca <- bca(median1.vec, conf.level = level)
                    median1.lb <- median1.bca[1]
                    median1.ub <- median1.bca[2]
                    exp.dur2.bca <- apply(exp.dur2.mat, 1, FUN=function(x){
                         bca(x, conf.level = level)
                    })
                    exp.dur2.lb <- exp.dur2.bca[1,]
                    exp.dur2.ub <- exp.dur2.bca[2,]
                    mean2.bca <- bca(mean2.vec, conf.level = level)
                    mean2.lb <- mean2.bca[1]
                    mean2.ub <- mean2.bca[2]
                    median2.bca <- bca(median2.vec, conf.level = level)
                    median2.lb <- median2.bca[1]
                    median2.ub <- median2.bca[2]
                    diff.bca <- apply(diff.mat, 1, FUN=function(x){
                         bca(x, conf.level = level)
                    })
                    diff.lb <- diff.bca[1,]
                    diff.ub <- diff.bca[2,]
                    mean.diff.bca <- bca(mean.diff.vec, conf.level = level)
                    mean.diff.lb <- mean.diff.bca[1]
                    mean.diff.ub <- mean.diff.bca[2]
                    median.diff.bca <- bca(median.diff.vec, conf.level = level)
                    median.diff.lb <- median.diff.bca[1]
                    median.diff.ub <- median.diff.bca[2]
               }
               res <- list(exp.dur1 = data.frame(exp.dur = exp.dur[,1], bootstrap.se = exp.dur1.se,
                                                 lb = exp.dur1.lb, ub = exp.dur1.ub),
                           mean1 = data.frame(mean = mean(exp.dur[,1], na.rm=TRUE), bootstrap.se = mean1.se,
                                              lb = mean1.lb, ub = mean1.ub),
                           median1 = data.frame(median = median(exp.dur[,1], na.rm=TRUE), bootstrap.se = median1.se,
                                                lb = median1.lb, ub = median1.ub),
                           exp.dur2 = data.frame(exp.dur = exp.dur[,2], bootstrap.se = exp.dur2.se,
                                                 lb = exp.dur2.lb, ub = exp.dur2.ub),
                           mean2 = data.frame(mean = mean(exp.dur[,2], na.rm=TRUE), bootstrap.se = mean2.se,
                                              lb = mean2.lb, ub = mean2.ub),
                           median2 = data.frame(median = median(exp.dur[,2], na.rm=TRUE), bootstrap.se = median2.se,
                                                lb = median2.lb, ub = median2.ub),
                           diff = data.frame(exp.dur = exp.dur[,3], bootstrap.se = diff.se,
                                             lb = diff.lb, ub = diff.ub),
                           mean.diff = data.frame(mean = mean(exp.dur[,3], na.rm=TRUE), bootstrap.se = mean.diff.se,
                                                  lb = mean.diff.lb, ub = mean.diff.ub),
                           median.diff = data.frame(median = median(exp.dur[,3], na.rm=TRUE), bootstrap.se = median.diff.se,
                                                    lb = median.diff.lb, ub = median.diff.ub))
          }

     }

     if(method=="gam"){
          res$gam.model <- gam.model
          res$gam.data <- gam.data
     }
     if(method=="npsf"){
          res$baseline.functions <- baseline.functions
     }
     if(!is.null(newdata2)){
          class(res) <- "coxedMargin"
     } else {
          class(res) <- "coxedExpdur"
     }
     return(res)
}
