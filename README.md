datalimited2: More stock assessment methods for data-limited fisheries
======================================================================

Installation
------------

The R package datalimited2 can then be installed from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("cfree14/datalimited2")
```

The package implements the following catch-only stock assessment models:

- cMSY from Froese et al. 2017; see `?cmsy2`
- Boosted regression tree approach (zBRT) from Zhou et al. 2017; see `?zbrt`
- Optimized catch-only model (OCOM) from Zhou et al. 2016; see `?ocom`
- Refined ORCS approach (rORCS) from Free et al. 2017; see `?rorcs`
- Superensemble approaches from Jensen & Free 2018; see `?super1`, `?super4`, `?super12`

Citation
--------

To cite the package, cite the authors of the original assessment method (shown below)
and the following:

Free CM (2018) datalimited2: More stock assessment methods for data-limited fisheries.
R package version 0.1.0. https://github.com/cfree14/datalimited2


References
----------

Free, C.M., O.P. Jensen, J. Wiedenmann, and J.J. Deroba. 2017. The refined ORCS approach: a catch-based method for estimating stock status and catch limits for data-poor fish stocks. Fisheries Research 193:60-70.

Froese, R., N. Demirel, G. Coro, K.M. Kleisner, and H. Winker. 2017. Estimating fisheries reference points from catch and resilience. Fish & Fisheries 18(3):506-526. 

Kleisner, K., D. Zeller, R. Froese, and D. Pauly. 2013. Using global catch data for inferences on the world’s marine fisheries. Fish & Fisheries 14(3):293-311.

Kleisner, K., and D. Pauly. 2011. Stock-status plots of fisheries for Regional Seas, in: Christensen, V., S. Lai, M.L.D. Palomares, D. Zeller, D. Pauly (Eds.), The State of Biodiversity and Fisheries in Regional Seas. Fisheries Centre Research Reports 19(3). University of British Columbia, Vancouver, Canada, pp. 37-40.

Zhou, S., A.E. Punt, Y. Ye, N. Ellis, C.M. Dichmont, M. Haddon, D.C. Smith, and A.D.M. Smith. 2017. Estimating stock depletion level from patterns of catch history. Fish & Fisheries 18(4):742-751.

Zhou, S., Z. Chen, C.M. Dichmont, N. Ellis, M. Haddon, A.E. Punt, A.D.M. Smith, D.C. Smith, and Y. Ye. 2016. An optimised catch-only assessment method for data poor fisheries. In: Catch-based methods for data-poor fisheries. Report to FAO. CSIRO, Brisbane, Australia.







