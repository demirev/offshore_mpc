library(tidyr)
library(Hmisc) # for weighted quantile
library(DescTools) # for Gini

source("R/functions/utils.R")


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
  age = "RA0300", #"RA0300_B", # RA0300 is not present for privacy reasons in some
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

Wealthfiles <- Dfiles %>% 
  seq_along %>%
  lapply(function(i){
    D <- Dfiles[[i]] %>%
      inner_join(Hfiles[[i]], by = c("surv_id","country")) %>%
      inner_join(Pfiles[[i]], # id in slighntly different format 
                 by = c("surv_id","country","reference_person")) %>% 
      replace_na_all %>%
      mutate(managed_accounts = as.numeric(managed_accounts)) %>% # convert type
      filter(age >= 25 & age <= 60) %>%
      #filter(age >= 6 & age < 13) %>% # 6 is 25-29, 13 is 60-64
      #filter(permanent_income_flag == 2 | country %in% c("FR","FI")) %>% # 2 is 'normal' income
      transmute(
        surv_id = surv_id,
        country = country,
        age = age,
        income = total_income,
        net_wealth = net_wealth,
        liquid_assets = deposits + mutual_funds + business_wealth +
          shares + managed_accounts + life_insurance,
        weight = household_weight,
        permanent_income_flag = permanent_income_flag
      )
  })
Wealthfiles_obs <- Wealthfiles[[1]]$country %>% table
Wealthfiles_summary <- Wealthfiles %>% lapply(summary)


# Create Summaries --------------------------------------------------------
# Gini by country --

# helper functions
calcGini <- function(dset, cntr = NULL, lower_cut = 0,
                     wealthvar = "net_wealth",
                     weightvar = "weight", cntryvar = "country") {
  
  # wrapper aroung Gini that does some filtering
  
  if (!is.null(cntr)) {
    dset <- dset[dset[cntryvar] == cntr, ]
  } 
  
  if (!is.null(lower_cut)) {
    dset <- dset[dset[wealthvar] >= lower_cut, ]
  }
  
  
  wealth <- dset[[wealthvar]]
  weight <- dset[[weightvar]]
  
  Gini(wealth, n = round(weight)) # point estimate
}


calcWealthPerc <- function(dset, cntr = NULL, probs = seq(0,1,0.05),
                           lower_cut = 0, wealthvar = "net_wealth",
                           weightvar = "weight", cntryvar = "country") {
  
  # weighted cumulative weight distribution quantiles

  if (!is.null(cntr)) {
    dset <- dset[dset[cntryvar] == cntr, ]
  } 
  
  if (!is.null(lower_cut)) {
    dset <- dset[dset[wealthvar] >= lower_cut, ]
  }
  
  dset <- dset[order(dset[[wealthvar]], decreasing = T), ] # descending wealth
  
  wealth <- dset[[wealthvar]]
  weight <- dset[[weightvar]]
  
  total_wealth <- sum(wealth)
  cumul_wealth <- cumsum(wealth)
  
  wtd.quantile(x = cumul_wealth, weights = weight, probs = probs) / total_wealth
  
}

calcWealthIncomeRatios <- function(
  dset, cntr = NULL, lower_cut = 0,
  probs = seq(0,1,0.05),
  wealthvar = "net_wealth",
  incomevar = "income",
  weightvar = "weight",
  cntryvar = "country") {
  
  # wealth to income ratios
  # EG: Poorest 10% have 2 w/i ratio, poorest 20% have ... etc
  
  if (!is.null(cntr)) {
    dset <- dset[dset[cntryvar] == cntr, ]
  } 
  
  if (!is.null(lower_cut)) {
    dset <- dset[dset[wealthvar] >= lower_cut, ]
  }
  
  dset <- dset[order(dset[[wealthvar]], decreasing = F), ] # !!!Ascending wealth
  
  wealth <- dset[[wealthvar]]
  weight <- dset[[weightvar]]
  income <- dset[[incomevar]]
  
  total_wealth <- sum(wealth)
  cumul_wealth <- cumsum(wealth)
  cumul_income <- cumsum(income)
  ratio  <- cumul_wealth/(cumul_income/4)
  
  qs <- wtd.quantile(x = cumul_wealth, weights = weight, probs = probs) / total_wealth
  
  sapply(probs, function(q) {
    ratio[cumul_wealth/total_wealth >= qs[which(probs==q)]][1] # first that meets 
    # the corresponding quantile of the wealth distribution
  })
  
}

# Ginis
Gini_netwealth <- sapply(unique(Wealthfiles[[1]]$country), function(cnt) {
  mi_point(Wealthfiles, FUN = function(dset) {
    calcGini(dset, cntr=cnt, wealthvar = "net_wealth")
  })
})

Gini_liquidassets <- sapply(unique(Wealthfiles[[1]]$country), function(cnt) {
  mi_point(Wealthfiles, FUN = function(dset) {
    calcGini(dset, cntr=cnt, wealthvar = "liquid_assets")
  })
})

# Wealth Quantiles
Prop_netwealth <- lapply(unique(Wealthfiles[[1]]$country), function(cnt) {
  mi_point(Wealthfiles, FUN = function(dset) {
    calcWealthPerc(dset, cntr=cnt, wealthvar = "net_wealth")
  })
})
names(Prop_netwealth) <- unique(Wealthfiles[[1]]$country)

Prop_liquidassets <- lapply(unique(Wealthfiles[[1]]$country), function(cnt) {
  mi_point(Wealthfiles, FUN = function(dset) {
    calcWealthPerc(dset, cntr=cnt, wealthvar = "liquid_assets")
  })
})
names(Prop_liquidassets) <- unique(Wealthfiles[[1]]$country)
 
# Wealth-to-income
WtoI_netwealth <- lapply(unique(Wealthfiles[[1]]$country), function(cnt) {
  mi_point(Wealthfiles, FUN = function(dset) {
    calcWealthIncomeRatios(dset, cntr=cnt, wealthvar = "net_wealth")
  })
})
names(WtoI_netwealth) <- unique(Wealthfiles[[1]]$country)

WtoI_liquidassets <- lapply(unique(Wealthfiles[[1]]$country), function(cnt) {
  mi_point(Wealthfiles, FUN = function(dset) {
    calcWealthIncomeRatios(dset, cntr=cnt, wealthvar = "liquid_assets")
  })
})
names(WtoI_liquidassets) <- unique(Wealthfiles[[1]]$country)


# Visualizations ----------------------------------------------------------


