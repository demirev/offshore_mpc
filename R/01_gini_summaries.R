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
      )
  })

# Create Summaries --------------------------------------------------------
allCountries <- unique(Wealthfiles[[1]]$country)

Gini_netwealth <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_Gini(
        dset,
        cntr=cnt,
        wealthvar = "net_wealth",
        filters = filter_both # age and non-negative
      )
    })
  })

Gini_liquidassets <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_Gini(
        dset,
        cntr=cnt,
        wealthvar = "liquid_assets",
        filters = filter_both # age and non-negative
      )
    })
  })

Prop_netwealth <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "net_wealth",
        filters = filter_both,  # age and non-negative
        descending = T # richest to poorest
      )
    })
  }) 
names(Prop_netwealth) <- allCountries

Prop_liquidassets <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "liquid_assets",
        filters = filter_both, # age and non-negative
        descending = T #richest to poorest
      )
    })
  })
names(Prop_liquidassets) <- allCountries


Prop_netwealth_asc <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "net_wealth",
        filters = filter_both, # age and non-negative
        descending = F # poorest to richerst
      )
    })
  }) 
names(Prop_netwealth_asc) <- allCountries

Prop_liquidassets_asc <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "liquid_assets",
        filters = filter_both, # age and non-negative
        descending = F #richest to poorest
      )
    })
  })
names(Prop_liquidassets_asc) <- allCountries

# Wealth-to-income
WtoI_netwealth <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WtoI(
        dset,
        cntr=cnt,
        wealthvar = "net_wealth",
        filters = filter_three,
        cumulative = F, # average per percentile
        probs = c(0,0.25,0.5,0.75,1),
        descending = F # poorest to richest
      )
    })
  })
names(WtoI_netwealth) <- allCountries

WtoI_liquidassets <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WtoI(
        dset,
        cntr=cnt,
        wealthvar = "liquid_assets",
        filters = filter_three,
        cumulative = F, # average per percentile
        probs = c(0,0.25,0.5,0.75,1),
        descending = F # poorest to richest
      )
    })
  })
names(WtoI_liquidassets) <- allCountries


# Visualizations ----------------------------------------------------------
v1 <- Prop_netwealth %>% as_tibble %>% lorenz_plotter
v2 <- tibble(
  net_wealth = Gini_netwealth,
  liquid_assets = Gini_liquidassets
) %>%
  bar_plotter(allCountries)
v3 <- WtoI_netwealth %>% box_plotter()


# dens <- with(Wealthfiles[[1]], tapply(net_wealth, INDEX = country, density))
# df <- data.frame(
#   x = unlist(lapply(dens, "[[", "x")),
#   y = unlist(lapply(dens, "[[", "y")),
#   country = rep(names(dens), each = length(dens[[1]]$x))
# )

# hist_plot <- plot_ly(df, x = ~x, y = ~y, color = ~country) %>%
#   add_lines() %>%
#   layout(
#     title = "Gini Coefficients",
#     barmode = 'group',
#     xaxis = list(
#       title = ""
#     ),
#     yaxis = list(
#       title = "",
#       showgrid = F
#     )
#   )


