#' Generate simulated durations using a baseline survivor function and proportional hazards
#'
#' This function is called by \code{\link[coxed]{sim.survdata}} and is not intended to be used by itself.
#' @param baseline The baseline hazard, cumulative hazard, survival, failure PDF, and failure CDF as output by \code{\link[coxed]{baseline.build}}
#' @param X A user-specified data frame containing the covariates that condition duration. If \code{NULL}, covariates are generated from
#' normal distributions with means given by the \code{mu} argument and standard deviations given by the \code{sd} argument
#' @param N Number of observations in each generated data frame
#' @param type If "none" (the default) data are generated with no time-varying covariates or coefficients.
#' If "tvc", data are generated with time-varying covariates, and if "tvbeta" data are generated with time-varying
#' coefficients (see details)
#' @param beta A user-specified vector containing the coefficients that for the linear part of the duration model. If \code{NULL}, coefficients are generated from
#' normal distributions with means of 0 and standard deviations of 0.1
#' @param xvars The number of covariates to generate. Ignored if \code{X} is not \code{NULL}
#' @param mu If scalar, all covariates are generated to have means equal to this scalar. If a vector, it specifies the mean of each covariate separately,
#' and it must be equal in length to \code{xvars}. Ignored if \code{X} is not \code{NULL}
#' @param sd If scalar, all covariates are generated to have standard deviations equal to this scalar. If a vector, it specifies the standard deviation
#' of each covariate separately, and it must be equal in length to \code{xvars}. Ignored if \code{X} is not \code{NULL}
#' @param censor The proportion of observations to designate as being right-censored
#' @details If \code{type="none"} then the function generates idiosyncratic survival functions for each observation via proportional hazards: first the
#' linear predictor is calculated from the X variables and beta coefficients, then the linear predictor is exponentiated and set as the exponent of the
#' baseline survivor function.  For each individual observation's survival function, a duration is drawn by drawing a single random number on U[0,1]
#' and finding the time point at which the survival function first decreases past this value. See Harden and Kropko (2018) for a more detailed description
#' of this algorithm.
#'
#' If \code{type="tvc"}, this function cannot accept user-supplied data for the covariates, as a time-varying covariate is expressed over time frames
#' which themselves convey part of the variation of the durations, and we are generating these durations. If user-supplied X data is provided, the
#' function passes a warning and generates random data instead as if \code{X=NULL}. Durations are drawn again using proportional hazards, and are passed
#' to the \code{\link[PermAlgo]{permalgorithm}} function in the \code{PermAlgo} package to generate the time-varying data structure (Sylvestre and Abrahamowicz 2008).
#'
#' If \code{type="tvbeta"} the first coefficient, whether coefficients are user-supplied or randomly generated, is interacted with the natural log of
#' the time counter from 1 to \code{T} (the maximum time point for the \code{baseline} functions). Durations are generated via proportional hazards,
#' and coefficients are saved as a matrix to illustrate their dependence on time.
#' @return Returns a list with the following components:
#' \tabular{ll}{
#' \code{data} \tab The simulated data frame, including the simulated durations, the censoring variable, and covariates\cr
#' \code{beta} \tab The coefficients, varying over time if \code{type} is "tvbeta" \cr
#' \code{XB} \tab The linear predictor for each observation \cr
#' \code{exp.XB} \tab The exponentiated linear predictor for each observation \cr
#' \code{survmat} \tab An (\code{N} x \code{T}) matrix containing the individual survivor function at
#' time t for the individual represented by row n   \cr
#' \code{tvc} \tab A logical value indicating whether or not the data includes time-varying covariates \cr
#' }
#' @references Harden, J. J. and Kropko, J. (2018). Simulating Duration Data for the Cox Model.
#' \emph{Political Science Research and Methods} \url{https://doi.org/10.1017/psrm.2018.19}
#'
#' Sylvestre M.-P., Abrahamowicz M. (2008) Comparison of algorithms to generate event times conditional on time-dependent covariates. \emph{Statistics in Medicine} \strong{27(14)}:2618–34.
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @seealso \code{\link[coxed]{sim.survdata}}, \code{\link[PermAlgo]{permalgorithm}}
#' @export
#' @examples baseline <- baseline.build(T=100, knots=8, spline=TRUE)
#' simdata <- generate.lm(baseline, N=1000, xvars=5, mu=0, sd=1, type="none", censor=.1)
#' summary(simdata$data)
#' simdata <- generate.lm(baseline, N=1000, xvars=5, mu=0, sd=1, type="tvc", censor=.1)
#' summary(simdata$data)
#' simdata <- generate.lm(baseline, N=1000, xvars=5, mu=0, sd=1, type="tvbeta", censor=.1)
#' simdata$beta
generate.lm <- function(baseline, X=NULL, N=1000, type="none", beta=NULL, xvars=3, mu=0, sd=1, censor=.1){
     T <- max(baseline$time)
     if(type=="none"){
          if(is.null(X)) X <- cbind(matrix(rnorm(N*xvars, mean=mu, sd=sd), N, xvars))
          if(!is.null(X)) X <- as.matrix(X)
          if(is.null(beta)) beta <- as.matrix(rnorm(ncol(X), mean=0, sd=.1))
          if(!is.null(beta)) beta <- as.matrix(beta)
          XB <- X%*%beta
          survival <- t(sapply(XB, FUN=function(x){baseline$survivor^exp(x)}, simplify=TRUE))
          survival <- cbind(1, survival)
          y <- apply(survival, 1, FUN=function(x){
               z <- diff(x < runif(1))
               r <- ifelse(all(z==0), T, which.max(z))
               return(r)
          })
          data <- data.frame(X)
          data$y <- y
          tvc <- FALSE
     } else if(type=="tvc"){
          if(!is.null(X)) warning("User-supplied X matrices are not implemented when type='tvc'. Generating random X data")
          if(is.null(beta)) beta <- as.matrix(rnorm(xvars, mean=0, sd=.1))
          if(!is.null(beta)) beta <- as.matrix(beta)
          Xmat1 <- expand.grid(N = 1:N, T=1:T)
          Xmat1$X1 <- rnorm(N*T)
          Xmat2 <- data.frame(matrix(rnorm(xvars*N), N, xvars))
          Xmat2$N <- 1:N
          Xmat2 <- Xmat2[,-1]
          X <- as.matrix(merge(Xmat1, Xmat2, by="N")[,-c(1,2)])
          XB <- matrix(X%*%beta, N, T, byrow=TRUE)
          survival <- t(apply(XB, 1, FUN=function(x){baseline$survivor^exp(x)}))
          survival <- cbind(1, survival)
          lifetimes <- apply(survival, 1, FUN=function(x){
               z <- diff(x < runif(1))
               r <- ifelse(all(z==0), T, which.max(z))
               return(r)
          })
          cen <- quantile(lifetimes, 1-censor) #Step 1: find the 1-censor quantile of the dist. of lifetimes
          m <- (2*cen + 1)/T #Step 2: set mean of the uniform distribution to this quantile. Anything above the uniform will be censored
          data <- PermAlgo::permalgorithm(N, T, X,
                                          XmatNames = colnames(X),
                                          eventRandom = lifetimes,
                                          censorRandom = round(runif(N, 1, m*T),0),
                                          betas = beta,
                                          groupByD = FALSE)
          data <- dplyr::rename(data, id=Id, failed=Event, start=Start, end=Stop)
          data <- dplyr::select(data, -Fup)
          rownames(data) <- NULL
          tvc <- TRUE
     } else if(type=="tvbeta"){
          if(is.null(X)) X <- matrix(rnorm(N*xvars, mean=mu, sd=sd), N, xvars)
          if(!is.null(X)) X <- as.matrix(X)
          if(is.data.frame(beta)) beta <- as.matrix(beta)
          if(!is.null(beta) & is.matrix(beta)){
                  d <- dim(beta)
                  if(d[1] != T) stop("User specified coefficients with type='tvbeta' must be specified either once, or for all T time points")
                  if(d[2] != xvars) stop("User specified coefficients must be specified for all X variables")
                  beta.mat <- beta
          }
          if(!is.null(beta) & !is.matrix(beta)){
                  if(length(beta)!=xvars) stop("There must be one user-supplied coefficient for each X variable")
                  beta.mat <- matrix(rep(beta, T), T, length(beta), byrow = TRUE)
                  beta.mat[,1] <- beta.mat[,1]*log(1:T)
          }
          if(is.null(beta)){
                  beta <- rnorm(ncol(X), mean=0, sd=.1)
                  beta.mat <- matrix(rep(beta, T), T, length(beta), byrow = TRUE)
                  beta.mat[,1] <- beta.mat[,1]*log(1:T)
          }

          XB <- apply(beta.mat, 1, FUN = function(b){
                  return(X %*% b)
          })
          survival <- t(apply(XB, 1, FUN=function(x){baseline$survivor^exp(x)}))
          survival <- cbind(1, survival)
          lifetimes <- apply(survival, 1, FUN=function(x){
               z <- diff(x < runif(1))
               r <- ifelse(all(z==0), T, which.max(z))
               return(r)
          })
          data <- data.frame(X)
          data <- dplyr::mutate(data, y = lifetimes)
          beta <- beta.mat
          tvc <- FALSE
     } else {stop("type must be one of 'none', 'tvc', or 'tvbeta'")}
     return(list(data=data, beta=beta, XB=XB, exp.XB = exp(XB), survmat=survival, tvc=tvc))
}
