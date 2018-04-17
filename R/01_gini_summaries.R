library(tidyr)
library(Hmisc) # for weighted quantile
library(DescTools) # for Gini

source("R/functions/utils.R")


# Prepare Data ------------------------------------------------------------
# import
Hfiles <- importType(dir = "data/HFCS_UDB_1_3_ASCII/", type = "H")
Dfiles <- importType(dir = "data/HFCS_UDB_1_3_ASCII/", type = "D")
Wfile  <- read_csv("data/HFCS_UDB_1_3_ASCII/W.csv")

var_dict <- list(
  id = "ID", # unqiue id across files
  surv_id = "SA0010", # id within survey
  country = "SA0100", # country id
  household_weight = "HW0010",
  net_wealth = "DN3001",
  deposits   = "DA2101",
  mutual_funds = "DA2102",
  business_wealth = "DA2104",
  shares = "DA2105",
  managed_accounts = "DA2106",
  life_insurance = "DA2109"
)

# only needed variables from each file
filter_dict <- function(dataset){
  dataset <- dataset[colnames(dataset) %in% var_dict]
  colnames(dataset) <- colnames(dataset) %>% 
    sapply(function(name) names(var_dict)[var_dict==name])
  return(dataset)
}

Hfiles <- Hfiles %>%
  lapply(filter_dict)

Dfiles <- Dfiles %>%
  lapply(filter_dict)

Wfile <- Wfile[colnames(Wfile) %in% 
                 c(var_dict, paste0("wr",str_extract(var_dict, "\\d+")))]
# note: Wfile contains no replication weights for derived variables (!)

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
      inner_join(Hfiles[[i]], by = c("id","surv_id","country")) %>%
      replace_na_all %>%
      mutate(managed_accounts = as.numeric(managed_accounts)) %>% # convert type
      transmute(
        id = id,
        surv_id = id,
        country = country,
        net_wealth = net_wealth,
        liquid_assets = deposits + mutual_funds + business_wealth +
          shares + managed_accounts + life_insurance,
        weight = household_weight
      )
  })

Wealthfiles_summary <- Wealthfiles %>% lapply(summary)


# Create Summaries --------------------------------------------------------
# Create wealth quantile table



