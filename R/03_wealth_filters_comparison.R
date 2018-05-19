#source("R/functions/descriptives.R")
library(tidyr)
library(DescTools) # for Gini
library(plotly)
library(kableExtra)

source("R/functions/utils.R")
source("R/functions/descriptives.R")

# Prepare Data ------------------------------------------------------------
# import
Wealthfiles = bigImport(targetdir = "data/HFCS_UDB_1_3_ASCII/")

# Calculate quintiles --------------------------------------------------------
allCountries <- unique(Wealthfiles[[1]]$country)

get_no_observations <- function(dset, wealthvar = "net_wealth",
                                filters = function(dset, ...){return(dset)}) {
  # returns a list with the number of observations for each country
  # in dset after applying filters using wealthvar
  dset <- filters(dset=dset, wealthvar=wealthvar)
  obs <- by(dset, dset$country, count)
  return(obs)
}


# this should really be in a loop!
quintiles <- list()

quintiles$raw <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        probs = seq(0, 1, 0.2) # quintiles
      )
    })
  })
quintiles$raw <- rbind(
  quintiles$raw,
  get_no_observations(Wealthfiles[[1]])
)

quintiles$age <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        probs = seq(0, 1, 0.2), # quintiles
        filters = filter_age # age filter only
      )
    })
  })
quintiles$age <- rbind(
  quintiles$age,
  get_no_observations(Wealthfiles[[1]], filters = filter_age)
)

quintiles$income <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        probs = seq(0, 1, 0.2), # quintiles
        filters = filter_income
      )
    })
  })
quintiles$income <- rbind(
  quintiles$income,
  get_no_observations(Wealthfiles[[1]], filters = filter_income)
)

quintiles$negative <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        probs = seq(0, 1, 0.2), # quintiles
        filters = filter_negative
      )
    })
  })
quintiles$negative <- rbind(
  quintiles$negative,
  get_no_observations(Wealthfiles[[1]], filters = filter_negative)
)

quintiles$age_negative <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        probs = seq(0, 1, 0.2), # quintiles
        filters = filter_both
      )
    })
  })
quintiles$age_negative <- rbind(
  quintiles$age_negative,
  get_no_observations(Wealthfiles[[1]], filters = filter_both)
)

quintiles$age_negative_income <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        probs = seq(0, 1, 0.2), # quintiles
        filters = filter_three
      )
    })
  })
quintiles$age_negative_income <- rbind(
  quintiles$age_negative_income,
  get_no_observations(Wealthfiles[[1]], filters = filter_three)
)

# save quintile data
#save(quintiles, file = "data/generated/quintiles_filters_12.RData")


# deciles for final calibration -------------------------------------------
allCountries <- unique(Wealthfiles[[1]]$country)

Targets <- list()

Targets$net <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "net_wealth",
        probs = seq(0, 1, 0.1), # deciles
        filters = filter_age, # age filter only
        negative_to_zero = F
      )
    })
  }) %>% t

Targets$net_offshore <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "offshore_wealth",
        probs = seq(0, 1, 0.1), # deciles
        filters = filter_age, # age filter only
        negative_to_zero = F
      )
    })
  }) %>% t

Targets$liq <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "liquid_assets",
        probs = seq(0, 1, 0.1), # deciles
        filters = filter_age, # age filter only
        negative_to_zero = F
      )
    })
  }) %>% t

Targets$liq_offshore <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "liquid_offshore_wealth",
        probs = seq(0, 1, 0.1), # deciles
        filters = filter_age, # age filter only
        negative_to_zero = F
      )
    })
  }) %>% t

# Output for calibration file
# This is very ugly
# dput(deciles_net, "")
# dput(deciles_net_off, "")
# dput(deciles_liq, "")
# dput(deciles_net_off, "")

# Save instead
save(Targets, file = 'data/generated/decile_targets.RData')