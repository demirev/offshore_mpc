#source("R/functions/descriptives.R")
library(tidyr)
library(DescTools) # for Gini
library(plotly)
library(kableExtra)

source("R/functions/utils.R")
source("R/functions/descriptives.R")

# Prepare Data ------------------------------------------------------------
# import
targetdir <- "data/HFCS_UDB_1_3_ASCII/"
Pfiles <- importType(dir = targetdir, type = "P")
Hfiles <- importType(dir = targetdir, type = "H")
Dfiles <- importType(dir = targetdir, type = "D")
Wfile  <- read_csv(paste0(targetdir,"W.csv"))

var_dict <- list(
  #id = "ID",
  surv_id = "SA0010", # id within survey
  country = "SA0100", # country id
  reference_person = "DHIDH1", # head of household id
  reference_person = "RA0010", # matches id of head of household from Pfiles
  household_weight = "HW0010",
  net_wealth = "DN3001",
  deposits   = "DA2101",
  mutual_funds = "DA2102",
  business_wealth = "DA2104",
  shares = "DA2105",
  managed_accounts = "DA2106",
  life_insurance = "DA2109",
  total_income = "DI2000",
  age = "RA0300",
  age_bin = "RA0300_B",# RA0300 is not present for privacy reasons in some
  permanent_income_flag = "HG0700"
)

# only needed variables from each file
filter_dict <- function(dataset){
  dataset <- dataset[colnames(dataset) %in% var_dict]
  colnames(dataset) <- colnames(dataset) %>% 
    sapply(function(name) names(var_dict)[var_dict==name])
  return(dataset)
}

Hfiles <- Hfiles %>%
  lapply(filter_dict) # for consistency

Dfiles <- Dfiles %>%
  lapply(filter_dict) %>%
  lapply(function(ds) ds[colnames(ds)!="household_weight"])

Pfiles <- Pfiles %>%
  lapply(filter_dict)

# summaries
Dfiles_summary <- Dfiles %>%
  lapply(summary)

# data prep - merge weights, NA to 0, derive total liquid assets
replace_na_all <- function(dataset) {
  replace_list <- lapply(dataset, function(x){return(0)})
  names(replace_list) <- colnames(dataset)
  replace_na(dataset, replace_list)
}

OffShoreWealth <- list(
  AT = 30.4679290417694,
  BE = 81.9924137825682,
  CY = 0,
  DE = 551.7731683678,
  ES = 163.79217913614,
  FI = 6.88646359030688,
  FR = 409.307786028679,
  GR = 115.166013805294,
  IT = 262.150676249415,
  LU = 0,
  MT = 0,
  NL = 50.6051185244171,
  PT = 50.9767183258964,
  SI = 1.34384376735097,
  SK = 2.09820099504479
) %>%
  lapply(function(x) x*0.7477) %>% # to eur
  lapply(function(x) x*1000000000) # from bln


Wealthfiles <- Dfiles %>% 
  seq_along %>%
  lapply(function(i){
    D <- Dfiles[[i]] %>%
      inner_join(Hfiles[[i]], by = c("surv_id","country")) %>%
      inner_join(Pfiles[[i]], by = c("surv_id","country","reference_person")) %>% 
      replace_na_all %>% # NA would become zeros here
      mutate(managed_accounts = as.numeric(managed_accounts)) %>% # convert type
      transmute(
        surv_id = surv_id,
        country = country,
        age = age,
        age_bin = age_bin,
        income = total_income,
        net_wealth = net_wealth,
        liquid_assets = deposits + mutual_funds + business_wealth +
          shares + managed_accounts + life_insurance,
        weight = household_weight,
        permanent_income_flag = permanent_income_flag
      ) %>%
      group_by(country) %>%
      mutate(
        offshore_wealth = calc_offshore2(
          wealthvar = net_wealth,#liquid_assets, 
          offshore = OffShoreWealth[country][[1]],
          weight = weight
        ),
        liquid_offshore_wealth = calc_offshore2(
          wealthvar = liquid_assets, 
          offshore = OffShoreWealth[country][[1]],
          weight = weight
        )
      ) %>%
      ungroup
  })

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
save(quintiles, file = "data/generated/quintiles_filters_13.RData")
