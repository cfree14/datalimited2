datalimited2: More stock assessment methods for data-limited fisheries
======================================================================

Installation
------------

Before installing datalimited2, you will need to install [JAGS](http://mcmc-jags.sourceforge.net) (for BSM) and a [C++ compiler](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) (i.e., Xcode for Macs, RTools for Windows; for installing packages using the "devtools" package).

The "datalimited2" R package can then be installed from GitHub with:

``` r
# Run if you don't already have devtools installed
install.packages("devtools")

# Run once devtools is successfully installed
devtools::install_github("cfree14/datalimited2")
library(datalimited2)
```

The package implements the following catch-only stock assessment models:

- cMSY from Froese et al. 2017; see `?cmsy2`
- Bayesian surplus production model (BSM) from Froese et al. 2017; see `?bsm`
- Boosted regression tree approach (zBRT) from Zhou et al. 2017a; see `?zbrt`
- Optimized catch-only model (OCOM) from Zhou et al. 2017b; see `?ocom`
- Refined ORCS approach (rORCS) from Free et al. 2017; see `?rorcs`
- Superensemble approaches from Jensen & Free 2018; see `?super1`, `?super4`, `?super12`

The results of all but the rORCS model can be visualized using the `?plot_dlm` function.

Citation
--------

To cite the package, cite the authors of the original assessment method (shown below)
and the following:

Free CM (2018) datalimited2: More stock assessment methods for data-limited fisheries.
R package version 0.1.0. https://github.com/cfree14/datalimited2


References
----------

Free CM, Jensen OP, Wiedenmann J, Deroba JJ (2017) The refined ORCS approach: a catch-based method for estimating stock status and catch limits for data-poor fish stocks. Fisheries Research 193: 60-70.

Froese R, Demirel N, Coro G, Kleisner KM, Winker H (2017) Estimating fisheries reference points from catch and resilience. Fish & Fisheries 18(3): 506-526. 

Zhou S, Punt AE, Ye Y, Ellis N, Dichmont CM, Haddon M, Smith DC, Smith ADM (2017a) Estimating stock depletion level from patterns of catch history. Fish & Fisheries 18(4): 742-751.

Zhou S, Punt AE, Smith ADM, Ye Y, Haddon M, Dichmont CM, Smith DC
(2017b) An optimised catch-only assessment method for data poor fisheries.
ICES Journal of Marine Science: doi:10.1093/icesjms/fsx226.







