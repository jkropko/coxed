#' Duration-Based Quantities of Interest and Simulation Methods for the Cox Proportional Hazards Model
#'
#' The Cox proportional hazards model (implemented in R with \code{\link[survival]{coxph}} in the \code{survival}
#' package or with \code{\link[rms]{cph}} in the \code{rms} package) is one of the most frequently used estimators
#' in duration (survival) analysis. Because it is estimated using only the observed durations' rank ordering, typical
#' quantities of interest used to communicate results of the Cox model come from the hazard function (e.g., hazard ratios
#' or percentage changes in the hazard rate). These quantities are substantively vague and difficult for many audiences
#' of research to understand. The \code{coxed} package introduces a suite of methods to address these problems.
#' The package allows researchers to calculate duration-based quantities from Cox model results, such as the expected
#' duration (or survival time) given covariate values and marginal changes in duration for a specified change in a
#' covariate. These duration-based quantities often match better with researchers' substantive interests and are
#' easily understood by most readers.
#' In addition, no standard method exists for simulating durations directly from the Cox model's data generating
#' process because it does not assume a distributional form for the baseline hazard function. The \code{coxed} package
#' also contains functions to simulate general duration data that does not rely on an assumption of any particular parametric
#' hazard function.
#' @section Duration-Based Quantities of Interest for the Cox Model:
#'
#' The \code{\link[coxed]{coxed}} function generates expected durations for individual
#' observations and/or marginal changes in expected duration given a change in a covariate
#' from the Cox proportional hazards model. Specifically, the methods can compute (1) the
#' expected duration for each observation used to fit the Cox model, given the covariates,
#' (2) the expected duration for a "new" observation with a covariate profile set by the
#' analyst, or (3) the first difference, or change, in expected duration given two new data frames.
#'
#' There are two different methods of generating duration-based quantities in the package.
#' \code{\link[coxed]{coxed}} with \code{type="npsf"} calculates expected durations by using the method proposed by
#' Cox and Oakes (1984, 107-109) for estimating the cumulative baseline hazard function.  This
#' method is nonparametric and results in a step-function representation of the cumulative
#' baseline hazard. Cox and Oakes (1984, 108) show that the cumulative baseline hazard function can be estimated
#' after fitting a Cox model by
#' \deqn{\hat{H}_0(t) = \sum_{\tau_j < t}\frac{d_j}{\sum_{l\in\Re(\tau_j)}\hat{\psi}(l)},}
#' where \eqn{\tau_j} represents time points earlier than \eqn{t}, \eqn{d_j} is a count of the
#' total number of failures at \eqn{\tau_j}, \eqn{\Re(\tau_j)} is the remaining risk set at \eqn{\tau_j},
#' and \eqn{\hat{\psi}(l)} represents the ELP from the Cox model for observations still in the
#' risk set at \eqn{\tau_j}. This equation is used calculate the cumulative baseline hazard at
#' all time points in the range of observed durations. This estimate is a stepwise function
#' because time points with no failures do not contribute to the cumulative hazard, so the function
#' is flat until the next time point with observed failures.
#'
#' We extend this method to obtain expected durations by first calculating the baseline survivor
#' function from the cumulative hazard function, using
#' \deqn{\hat{S}_0(t) = \exp[-\hat{H}_0(t)].}
#' Each observation's survivor function is related to the baseline survivor function by
#' \deqn{\hat{S}_i(t) = \hat{S}_0(t)^{\hat{\psi}(i)},}
#' where \eqn{\hat{\psi}(i)} is the exponentiated linear predictor (ELP) for observation \eqn{i}.
#' These survivor functions can be used directly to calculate expected durations for each
#' observation.  The expected value of a non-negative random variable can be calculated by
#' \deqn{E(X) = \int_0^{\infty} \bigg(1 - F(t)\bigg)dt,}
#' where \eqn{F(.)} is the cumulative distribution function for \eqn{X}.  In the case of a
#' duration variable \eqn{t_i}, the expected duration is
#' \deqn{E(t_i) = \int_0^T S_i(t)\,dt,}
#' where \eqn{T} is the largest possible duration and \eqn{S(t)} is the individual's survivor
#' function.  We approximate this integral with a right Riemann-sum by calculating the survivor
#' functions at every discrete time point from the minimum to the maximum observed durations,
#' and multiplying these values by the length of the interval between time points with observed failures:
#' \deqn{E(t_i) \approx \sum_{t_j \in [0,T]} (t_j - t_{j-1})S_i(t_j).}
#'
#' \code{\link[coxed]{coxed}} with \code{type="gam"} employs a generalized additive model (GAM) to map the model's estimated linear
#' predictor values to duration times and proceeds according to five steps. First, it uses coefficient
#' estimates from the Cox model, so researchers must first estimate the model just as they always have.
#' Then the method computes expected values of risk for each observation by matrix-multiplying the
#' covariates by the estimated coefficients from the model, then exponentiating the result. This creates
#' the exponentiated linear predictor (ELP). Then the observations are ranked from smallest to largest
#' according to their values of the ELP. This ranking is interpreted as the expected order of failure;
#' the larger the value of the ELP, the sooner the model expects that observation to fail, relative to
#' the other observations. The next step is to connect the model's expected risk for each observation (ELP) to duration time
#' (the observed durations). A \code{\link[mgcv]{gam}} fits a model to data by using a series of locally-estimated polynomial
#' splines set by the user (see, for example, Wood, Pya, and Saefken 2016). It is a flexible means of allowing for
#' the possibility of nonlinear relationships between variables. \code{coxed} with \code{type="gam"} uses a GAM to model the observed
#' durations as a function of the linear predictor ranks generated in the previous step. More specifically, the method
#' utilizes a cubic regression spline to draw a smoothed line summarizing the bivariate relationship between
#' the observed durations and the ranks. The GAM fit can be used directly to compute expected durations, given the covariates, for each observation
#' in the data.
#'
#' See Kropko and Harden (2018) for further details about generating expected durations and marginal changes in expected
#' duration from the Cox model. The \code{\link[coxed]{coxed}} function can also generate these quantities from data with time-varying
#' covariates (see \code{\link[coxed]{coxed.npsf.tvc}} and \code{\link[coxed]{coxed.gam.tvc}}).
#'
#' @section Simulating duration data for the Cox model:
#' The \code{\link[coxed]{sim.survdata}} function generates simulated duration data. It can accept a user-supplied
#' hazard function, or else it uses the flexible-hazard method described in Harden and Kropko (2018) to generate
#' a hazard that does not necessarily conform to any parametric hazard function. It can generate data with time-varying
#' covariates or coefficients. For time-varying covariates \code{type="tvc"} it employs the permutational algorithm by Sylvestre and Abrahamowicz (2008).
#' For time-varying coefficients with \code{type="tvbeta"}, the first beta coefficient that is either supplied by the user or generated by
#' the function is multiplied by the natural log of the failure time under consideration.
#'
#' The flexible-hazard method employed when \code{hazard.fun} is \code{NULL} generates a unique baseline hazard by fitting a curve to
#' randomly-drawn points. This produces a wide variety
#' of shapes for the baseline hazard, including those that are unimodal, multimodal, monotonically increasing or decreasing, and many other
#' shapes. The method then generates a density function based on each baseline hazard and draws durations from it in a way that circumvents
#' the need to calculate the inverse cumulative baseline hazard. Because the shape of the baseline hazard can vary considerably, this approach
#' matches the Cox model’s inherent flexibility and better corresponds to the assumed data generating process (DGP) of the Cox model. Moreover,
#' repeating this process over many iterations in a simulation produces simulated samples of data that better reflect the considerable
#' heterogeneity in data used by applied researchers. This increases the generalizability of the simulation results. See Harden and Kropko (2018)
#' for more detail.
#'
#' When generating a marginal effect, first the user specifies a covariate by typing its column number in the \code{X} matrix into the \code{covariate}
#' argument, then specifies the high and low values at which to fix this covariate.  The function calculates the differences in expected duration for each
#' observation when fixing the covariate to the high and low values.  If \code{compare} is \code{median}, the function reports the median of these differences,
#' and if \code{compare} is \code{mean}, the function reports the median of these differences, but any function may be employed that takes a vector as input and
#' outputs a scalar.
#'
#' If \code{censor.cond} is \code{FALSE} then a proportion of the observations specified by \code{censor} is randomly and uniformly selected to be right-censored.
#' If \code{censor.cond} is \code{TRUE} then censoring depends on the covariates as follows: new coefficients are drawn from normal distributions with mean 0 and
#' standard deviation of 0.1, and these new coefficients are used to create a new linear predictor using the \code{X} matrix.  The observations with the largest
#' (100 x \code{censor}) percent of the linear predictors are designated as right-censored.
#'
#' Finally, \code{link[coxed]{survsim.plot}} is useful for visualizing the baseline functions, including hazard, that result from
#' \code{link[coxed]{sim.survdata}} for a particular draw of simulated duration data. The function uses \code{\link[ggplot2]{ggplot}}
#' to create line plots for the baseline failure PDF, the baseline failure CDF, the baseline survivor function,
#' and the baseline hazard function over time.  The baseline functions and time are attributes of the
#' \code{\link[coxed]{sim.survdata}} output.
#' @import survival
#' @import rms
#' @import mgcv
#' @importFrom stats aggregate end glm.fit median naresid predict pt pnorm qnorm quantile rnorm runif sd start time ts
#' @importFrom utils globalVariables
#' @examples
#' ## See the examples for coxed, sim.survdata, and survsim.plot
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden2@@nd.edu>
#' @references
#' Harden, J. J. and Kropko, J. (2018) Simulating Duration Data for the Cox Model.
#' \emph{Political Science Research and Methods} \url{https://doi.org/10.1017/psrm.2018.19}
#'
#' Kropko, J. and Harden, J. J. (2018) Beyond the Hazard Ratio: Generating Expected
#' Durations from the Cox Proportional Hazards Model. \emph{British Journal of Political Science}
#' \url{https://doi.org/10.1017/S000712341700045X}
#'
#' Box-Steffensmeier, J. M. (1996)
#' A Dynamic Analysis of The Role of War Chests in Campaign Strategy.
#' \emph{American Journal of Political Science} \strong{40} 352-371
#' \url{https://doi.org/10.2307/2111628}
#'
#' Hyman, J. M. (1983) Accurate monotonicity preserving cubic interpolation. \emph{SIAM J. Sci. Stat. Comput.} \strong{4}, 645–654.
#' \url{https://doi.org/10.1137/0904045}
#'
#' Martin, L. W and Vanberg, G. (2003) Wasting Time?
#' The Impact of Ideology and Size on Delay in Coalition Formation.
#' \emph{British Journal of Political Science} \strong{33} 323-344
#' \url{https://doi.org/10.1017/S0007123403000140}
#'
#' Sylvestre M.-P., Abrahamowicz M. (2008) Comparison of algorithms to generate event times conditional on time-dependent
#' covariates. \emph{Statistics in Medicine} \strong{27(14)}:2618–34. \url{https://doi.org/10.1002/sim.3092}
#'
#' Wood, S.N., N. Pya and B. Saefken (2016) Smoothing parameter and model selection for general smooth models (with discussion). \emph{Journal of the American Statistical Association} \strong{111}, 1548-1575
#' \url{http://dx.doi.org/10.1080/01621459.2016.1180986}
#' @docType package
#' @name coxed-package
NULL
