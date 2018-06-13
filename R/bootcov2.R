#' Boostrapping algorithm for \code{coxed}
#'
#' This function uses bootstrapping to create standard errors and
#' confidence intervals for the quantities produced by the coxed()
#' function. It is adapted from the \code{\link[rms]{bootcov}} function
#' in the \code{\link[rms]{rms}} package. It is called by the \code{\link[coxed]{coxed}}
#' function and is not intended to be used by itself.  Please refer to
#' the original \code{\link[rms]{bootcov}} function for general bootstrapping
#' applications.
#' @param fit an estimated Cox proportional hazards model object with
#' class "coxph" or "cph"
#' @param cluster a variable indicating groupings. \code{cluster} may
#' be any type of vector (factor, character, integer). Unique values
#' of \code{cluster} indicate possibly correlated groupings of observations.
#' Note the data used in the fit and stored in \code{fit$x} and \code{fit$y}
#' may have had observations containing missing values deleted. It is assumed
#' that if there were any \code{NA}s, an \code{naresid} function exists for
#' the class of fit. This function restores \code{NA}s so that the rows of
#' the design matrix coincide with \code{cluster}
#' @param B Number of bootstrap simulation iterations
#' @param fitter the name of a function with arguments \code{(x,y)} that will fit
#' bootstrap samples. Default is taken from the class of fit if it is \code{ols},
#' \code{lrm}, \code{cph}, \code{psm}, \code{Rq}. If \code{fitter="tvc"} the
#' function employs \code{\link[survival]{agreg.fit}}
#' @param coef.reps set to \code{TRUE} if you want to store a matrix of all
#' bootstrap regression coefficient estimates in the returned component \code{boot.Coef}.
#' @param loglik set to \code{TRUE} to store -2 log likelihoods for each bootstrap model,
#' evaluated against the original \code{x} and \code{y} data. The default is to do this when
#' \code{coef.reps} is specified as \code{TRUE}. The use of \code{loglik=TRUE} assumes
#' that an \code{oos.loglik} method exists for the type of model being analyzed, to calculate
#' out-of-sample -2 log likelihoods (see \code{rmsMisc}). After the B -2 log likelihoods
#' (stored in the element named \code{boot.loglik} in the returned fit object), the \code{B+1}
#' element is the -2 log likelihood for the original model fit
#' @param pr set to \code{TRUE} to print the current sample number to monitor progress
#' @param maxit maximum number of iterations, to pass to \code{fitter}
#' @param group a grouping variable used to stratify the sample upon bootstrapping. This
#' allows one to handle \code{k}-sample problems, i.e., each bootstrap sample will be forced
#' to select the same number of observations from each level of group as the number appearing
#' in the original dataset. You may specify both \code{group} and \code{cluster}
#' @param stat a single character string specifying the name of a \code{stats} element produced
#' by the fitting function to save over the bootstrap repetitions. The vector of saved statistics
#' will be in the \code{boot.stats} part of the list returned by \code{bootcov}
#'
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden@@nd.edu>, based
#' on the code for the \code{\link[rms]{bootcov}} function in the \code{\link[rms]{rms}} package
#' by Frank Harrell and Bill Pikounis
#'
#' @details This function contains the same code as the \code{\link[rms]{bootcov}} function in
#' the \code{\link[rms]{rms}} package, with a few alterations to work better with the \code{\link[coxed]{coxed}}
#' function. First, we output a result attribute \code{b.ind}, which contains the observation numbers from the estimation sample
#' that are drawn with replacement to produce the bootstrap sample and takes into account clustering.
#' Second, we program a new class, \code{tvc}, for
#' \code{fitter} to use \code{\link[survival]{agreg.fit}} instead of \code{\link[survival]{coxph.fit}}
#' when the data contain time-varying covariates.
#'
#' @seealso \code{\link[coxed]{coxed}}, \code{\link[survival]{coxph}}, \code{\link[rms]{cph}}, \code{\link[rms]{bootcov}}
#' @export
bootcov2 <- function (fit, cluster, B = 200, fitter, coef.reps = TRUE, loglik = FALSE,
    pr = FALSE, maxit = 15, group = NULL, stat = NULL)
{
    if(missing(cluster)) b.ind <- matrix(NA, nrow = nrow(fit$x), ncol = B)
    else b.ind <- list()
    coxcph <- inherits(fit, "coxph") || inherits(fit, "cph") ||
        (length(fit$fitFunction) && any(c("cph", "coxph") %in%
            fit$fitFunction))
    nfit <- fit$fitFunction[1]
    if (!length(nfit))
        nfit <- setdiff(oldClass(fit), "Design")[1]
    if (length(fit$weights) && (coxcph || nfit[1] == "Rq"))
        stop("does not handle weights")
    if (!length(X <- fit$x) | !length(Y <- fit$y))
        stop("you did not specify x=TRUE and y=TRUE in the fit")
    sc.pres <- match("scale", names(fit), 0) > 0
    ns <- fit$non.slopes
    if (nfit == "psm") {
        fixed <- fit$fixed
        fixed <- if (length(fixed) == 1 && is.logical(fixed) &&
            !fixed)
            list()
        else list(scale = TRUE)
        fixed <- NULL
        dist <- fit$dist
        parms <- fit$parms
    }
    if (nfit == "Glm")
        fitFamily <- fit$family
    penalty.matrix <- fit$penalty.matrix
    if (missing(fitter)) {
        fitter <- switch(nfit, ols = if (length(penalty.matrix)) {
            function(x, y, penalty.matrix, ...) {
                lm.pfit(x, y, penalty.matrix = penalty.matrix,
                  tol = 1e-11, regcoef.only = TRUE)
            }
        } else function(x, y, ...) {
            lm.fit.qr.bare(x, y, tolerance = 1e-11, intercept = FALSE)
        }, lrm = function(x, y, maxit = 15, penalty.matrix, ...) {
            lrm.fit(x, y, maxit = maxit, tol = 1e-11, penalty.matrix = penalty.matrix)
        }, cph = function(x, y, strata = NULL, maxit = 15, ...) {
            coxph.fit(x, y,
                      strata = strata,
                      method = "efron",
                      rownames = rownames(x),
                      control = coxph.control(eps = 1e-04,
                                        toler.chol = 1e-11,
                                        toler.inf=1,
                                        iter.max = maxit))
        }, tvc = function(x, y, strata = NULL, init=NULL, maxit = 15, ...) {
            agreg.fit(x, y,
                      strata = strata,
                      init = init,
                      method = "efron",
                      rownames = rownames(x),
                      control = coxph.control(eps = 1e-04,
                                        toler.chol = 1e-11,
                                        toler.inf=1,
                                        iter.max = maxit))
        }, psm = function(x, y, maxit = 15, ...) {
            survreg.fit2(x, y, dist = dist, parms = parms, fixed = fixed,
                offset = NULL, init = NULL, maxiter = maxit)
        }, bj = function(x, y, maxit = 15, eps = 1e-04, ...) {
            bj.fit(x, y, control = list(iter.max = maxit, eps = 1e-04))
        }, Glm = function(x, y, ...) {
            glm.fit(x, as.vector(y), family = fitFamily)
        }, Rq = RqFit(fit, wallow = FALSE))
    }
    if (!length(fitter))
        stop("fitter not valid")
    if (loglik) {
        oosl <- switch(nfit, ols = oos.loglik.ols, lrm = oos.loglik.lrm,
            cph = oos.loglik.cph, psm = oos.loglik.psm, Glm = oos.loglik.Glm)
        if (!length(oosl))
            stop("loglik=TRUE but no oos.loglik method for model in rmsMisc")
        Loglik <- double(B + 1)
        Loglik[B + 1] <- oosl(fit)
    }
    else Loglik <- NULL
    n <- nrow(X)
    p <- length(fit$coef)
    vname <- names(fit$coef)
    if (sc.pres) {
        p <- p + 1
        vname <- c(vname, "log scale")
    }
    fillInTheBlanks <- function(S) {
        L <- !is.na(S)
        c(S[L][1], S[L])[cumsum(L) + 1]
    }
    fill <- function(cof, vn, ns) {
        p <- length(vn)
        if (length(cof) == p)
            return(cof)
        nc <- ns - (p - length(cof))
        cints <- cof[1:nc]
        ints <- rep(NA, ns)
        names(ints) <- vn[1:ns]
        ints[names(cints)] <- cints
        if (is.na(ints[ns])) {
            l <- ns
            if (ns > 1) {
                for (j in (ns - 1):1) {
                  if (!is.na(ints[j]))
                    break
                  l <- j
                }
            }
            ints[l:ns] <- -Inf
        }
        c(rev(fillInTheBlanks(rev(ints))), cof[-(1:nc)])
    }
    bar <- rep(0, p)
    cov <- matrix(0, nrow = p, ncol = p, dimnames = list(vname,
        vname))
    if (coef.reps)
        coefs <- matrix(NA, nrow = B, ncol = p, dimnames = list(NULL,
            vname))
    if (length(stat))
        stats <- numeric(B)
    #Y <- as.matrix(if (is.category(Y))
    #    oldUnclass(Y)
    #else Y)
    ny <- ncol(Y)
    Strata <- fit$Strata
    nac <- fit$na.action
    if (length(group)) {
        if (length(group) > n) {
            if (length(nac)) {
                j <- !is.na(naresid(nac, Y) %*% rep(1, ny))
                group <- group[j]
            }
        }
        if (length(group) != n)
            stop("length of group does not match # rows used in fit")
        group.inds <- split(1:n, group)
        ngroup <- length(group.inds)
    }
    else ngroup <- 0
    anyinf <- FALSE
    if (missing(cluster)) {
        b <- 0
        for (i in 1:B) {
            if (pr)
                cat(i, "\r")
            if (ngroup) {
                j <- integer(n)
                for (si in 1:ngroup) {
                  gi <- group.inds[[si]]
                  j[gi] <- sample(gi, length(gi), replace = TRUE)
                }
            }
            else j <- sample(1:n, n, replace = TRUE); b.ind[ , i] <- j
            f <- tryCatch(fitter(X[j, , drop = FALSE], Y[j, ,
                drop = FALSE], maxit = maxit, penalty.matrix = penalty.matrix,
                strata = Strata[j]), error = function(...) list(fail = TRUE))
            if (length(f$fail) && f$fail)
                next
            cof <- f$coefficients
            if (any(is.na(cof)))
                next
            b <- b + 1
            if (sc.pres)
                cof <- c(cof, `log scale` = log(f$scale))
            cof <- fill(cof, vname, ns)
            if (any(is.infinite(cof)))
                anyinf <- TRUE
            if (coef.reps)
                coefs[b, ] <- cof
            if (length(stat))
                stats[b] <- f$stats[stat]
            bar <- bar + cof
            cof <- as.matrix(cof)
            cov <- cov + cof %*% t(cof)
            if (loglik)
                Loglik[b] <- oosl(f, matxv(X, cof), Y)
        }
        if (pr)
            cat("\n")
    }
    else {
        if (length(cluster) > n) {
            if (length(nac)) {
                j <- !is.na(naresid(nac, Y) %*% rep(1, ny))
                cluster <- cluster[j]
            }
        }
        if (length(cluster) != n)
            stop("length of cluster does not match # rows used in fit")
        if (any(is.na(cluster)))
            stop("cluster contains NAs")
        cluster <- as.character(cluster)
        clusters <- unique(cluster)
        nc <- length(clusters)
        Obsno <- split(1:n, cluster)
        b <- 0
        for (i in 1:B) {
            if (pr)
                cat(i, "\r")
            if (ngroup) {
                j <- integer(0)
                for (si in 1:ngroup) {
                  gi <- group.inds[[si]]
                  cluster.gi <- cluster[gi]
                  clusters.gi <- unique(cluster.gi)
                  nc.gi <- length(clusters.gi)
                  Obsno.gci <- split(gi, cluster.gi)
                  j.gci <- sample(clusters.gi, nc.gi, replace = TRUE)
                  obs.gci <- unlist(Obsno.gci[j.gci])
                  j <- c(j, obs.gci)
                }
                obs <- j
            }
            else {
                j <- sample(clusters, nc, replace = TRUE)
                obs <- unlist(Obsno[j])
                b.ind[[i]] <- obs
            }
            f <- tryCatch(fitter(X[obs, , drop = FALSE], Y[obs,
                , drop = FALSE], maxit = maxit, penalty.matrix = penalty.matrix,
                strata = Strata[obs]), error = function(...) list(fail = TRUE))
            if (length(f$fail) && f$fail)
                next
            cof <- f$coefficients
            if (any(is.na(cof)))
                next
            b <- b + 1
            if (sc.pres)
                cof <- c(cof, `log scale` = log(f$scale))
            cof <- fill(cof, vname, ns)
            if (any(is.infinite(cof)))
                anyinf <- TRUE
            if (coef.reps)
                coefs[b, ] <- cof
            if (length(stat))
                stats[b] <- f$stats[stat]
            bar <- bar + cof
            cof <- as.matrix(cof)
            cov <- cov + cof %*% t(cof)
            if (loglik)
                Loglik[b] <- oosl(f, matxv(X, cof), Y)
        }
        if (pr)
            cat("\n")
    }
    if (b < B) {
        warning(paste("fit failure in", B - b, "resamples.  Might try increasing maxit"))
        if (coef.reps)
            coefs <- coefs[1:b, , drop = FALSE]
        Loglik <- Loglik[1:b]
    }
    if (anyinf)
        warning("at least one resample excluded highest Y values, invalidating bootstrap covariance matrix estimate")
    bar <- bar/b
    fit$B <- b
    names(bar) <- vname
    fit$boot.coef <- bar
    if (coef.reps)
        fit$boot.Coef <- coefs
    bar <- as.matrix(bar)
    cov <- (cov - b * bar %*% t(bar))/(b - 1)
    fit$orig.var <- fit$var
    fit$var <- cov
    fit$boot.loglik <- Loglik
    if(missing(cluster)) fit$b.ind <- b.ind
    else fit$b.ind <- as.matrix(unname(do.call(cbind, lapply(b.ind, ts))))
    if (length(stat))
        fit$boot.stats <- stats
    if (nfit == "Rq") {
        newse <- sqrt(diag(cov))
        newt <- fit$summary[, 1]/newse
        newp <- 2 * (1 - pt(abs(newt), fit$stats["n"] - fit$stats["p"]))
        fit$summary[, 2:4] <- cbind(newse, newt, newp)
    }
    fit
}

