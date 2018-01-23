
# Packages
library(plyr)
library(dplyr)

# User data
################################################################################

# cMSY example data
cmsy_data <- read.csv("data-raw/O_ICES&BS_Catch_7.csv", as.is=T)
SOLIRIS <- cmsy_data %>%
  filter(Stock=="sol-iris")

# OCOM example data
ocom_data <- read.csv("data-raw/ocom_data.csv", as.is=T)
TIGERFLAT <- ocom_data %>%
  filter(stock=="Tiger Flathead")

# Example RAMLDB stock
load("data-raw/YELLSNEMATL.Rdata")
YELLSNEMATL_orig <- YELLSNEMATL
YELLSNEMATL <- YELLSNEMATL_orig %>%
  select(stockid, year, tc, ssb, f, b_bmsy, f_fmsy) %>%
  rename(catch=tc, biomass=ssb, bbmsy=b_bmsy, ffmsy=f_fmsy)

# Build data
devtools::use_data(SOLIRIS, overwrite=T)
devtools::use_data(TIGERFLAT, overwrite=T)
devtools::use_data(YELLSNEMATL, overwrite=T)

# Internal data
################################################################################

# Load data
load("data-raw/BRTmodelP8.Rdata")
load("data-raw/BRTmodelP38.Rdata")
load("data-raw/refined_orcs_approach_bct_model.Rdata")

# Build data
devtools::use_data(BRTmodelP8, BRTmodelP38, rorcs_model, overwrite=T, internal=T)








