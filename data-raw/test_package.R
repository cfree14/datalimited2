
# Install package
devtools::install_github("cfree14/datalimited2")
library(datalimited2)

# Test BSM
output <- bsm(year=SOLIRIS$yr, catch=SOLIRIS$ct, biomass=SOLIRIS$bt, btype="CPUE", r.low=0.18, r.hi=1.02)
plot_dlm(output)

# Test rORCS
scores <- c(1, 2, NA, 2, 2, 3, 1.93, 2, 1, 2, 1, 3)
rorcs(scores)

# Test zBRT
output <- zbrt(year=TIGERFLAT$yr, catch=TIGERFLAT$catch)
plot_dlm(output)

# Test OCOM
output <- ocom(year=TIGERFLAT$yr, catch=TIGERFLAT$catch, m=0.27)
plot_dlm(output)

# Test cMSY
output <- cmsy2(year=SOLIRIS$yr, catch=SOLIRIS$ct, r.low=0.18, r.hi=1.02)
plot_dlm(output)







