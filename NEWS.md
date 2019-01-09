<!-- README.md is generated from README.Rmd. Please edit that file -->
coxed 0.2.3
===========
-  fixed bug with user-supplied hazards. If the hazard does not equal 0 at the latest time point, some number of observations can be right-censored naturally. Previously, the bug set these observations to the earliest, not latest time point. The bug is corrected and a warning is supplied

coxed 0.2.2
===========

-   fixed a bug with the TVC routines. All observations were being used to calculate cumulative baseline hazard for the NPSF function, resulting in a larger risk set, leading to lower hazard and longer durations than are correct. We fixed the bug by calculating the CBH from only the last observation of each unique value of ID

coxed 0.2.1
===========

-   edited and rebuilt vignettes 8/24/2018

coxed 0.2.0
===========

-   updated version number to reflect the added feature, and fixed minor bug 8/23/2018

coxed 0.1.3
===========

-   adds bias-corrected and accelerated (DiCiccio and Efron 1996) confidence intervals for bootstrapped samples in addition to the studentized and empirical options 8/23/2018

coxed 0.1.2
===========

-   changes a sentence in the coxed.Rmd vignette, and updated the version number for the next update 6/15/2018 \# coxed 0.1.1

-   first working version, ready for extensions/debugging/meeting CRAN demands 6/13/2018

-   changes for CRAN resubmission: reduced the coxed() example to run in half the time, and reduced the computational requirements of the vignettes 6/14/2018

coxed 0.1.0
===========

-   initial release 6/13/2018
