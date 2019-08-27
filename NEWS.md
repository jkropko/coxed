<!-- README.md is generated from README.Rmd. Please edit that file -->
coxed 0.3.1
===========

-   fixed a bug that led to (very close to 0) negative values of the
    survivor function when the random spline failure CDF reaches 1 prior
    to the final time point

coxed 0.3.0
===========

-   edits to help documentation and vignettes in preparation for release
    to CRAN

coxed 0.2.7
===========

-   fixed bug with coxed.npsf.tvc() in which all prior time points for
    every observation were included in the risk set. Changed it so that
    all yet-to-fail observations are included using ONLY the time point
    in question

coxed 0.2.6
===========

-   rewrote rank.predict() to make GAM much faster with out of sample
    prediction with GAM

coxed 0.2.5
===========

-   fixed a bug in cox.gam.tvc.R and cox.npsf.tvc.R in which the
    functions failed to produce out of sample predictions

coxed 0.2.4
===========

-   fixed typos in vignettes

coxed 0.2.3
===========

-   fixed a bug with sim.survdata() with type=“tvc” that prevented
    properly setting the amount of right-censoring

coxed 0.2.2
===========

-   fixed a bug with the TVC routines. All observations were being used
    to calculate cumulative baseline hazard for the NPSF function,
    resulting in a larger risk set, leading to lower hazard and longer
    durations than are correct. We fixed the bug by calculating the CBH
    from only the last observation of each unique value of ID

coxed 0.2.1
===========

-   edited and rebuilt vignettes 8/24/2018

coxed 0.2.0
===========

-   updated version number to reflect the added feature, and fixed minor
    bug 8/23/2018

coxed 0.1.3
===========

-   adds bias-corrected and accelerated (DiCiccio and Efron 1996)
    confidence intervals for bootstrapped samples in addition to the
    studentized and empirical options 8/23/2018

coxed 0.1.2
===========

-   changes a sentence in the coxed.Rmd vignette, and updated the
    version number for the next update 6/15/2018 \# coxed 0.1.1

-   first working version, ready for extensions/debugging/meeting CRAN
    demands 6/13/2018

-   changes for CRAN resubmission: reduced the coxed() example to run in
    half the time, and reduced the computational requirements of the
    vignettes 6/14/2018

coxed 0.1.0
===========

-   initial release 6/13/2018
