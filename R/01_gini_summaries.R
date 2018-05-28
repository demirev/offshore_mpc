#source("R/functions/descriptives.R")
library(tidyr)
library(DescTools) # for Gini
library(plotly)
library(kableExtra)

source("R/functions/utils.R")
source("R/functions/descriptives.R")

# Prepare Data ------------------------------------------------------------
# import
Wealthfiles <- bigImport()

# Create Summaries --------------------------------------------------------
allCountries <- unique(Wealthfiles[[1]]$country)

Totals <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_Totals(
        dset,
        cntr=cnt,
        wealthvar = "net_wealth", # for filter
        filters = filter_both # age and non-negative
      )
    })
  })

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

Gini_offshore <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_Gini(
        dset,
        cntr=cnt,
        wealthvar = "offshore_wealth",
        filters = filter_both # age and non-negative
      )
    })
  })

Gini_offshore_liquid <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_Gini(
        dset,
        cntr=cnt,
        wealthvar = "liquid_offshore_wealth",
        filters = filter_both # age and non-negative
      )
    })
  })

Gini_offshore_pareto <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_Gini(
        dset,
        cntr=cnt,
        wealthvar = "offshore_wealth_pareto",
        filters = filter_both # age and non-negative
      )
    })
  })

Gini_offshore_liquid_pareto <- allCountries %>%
  sapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_Gini(
        dset,
        cntr=cnt,
        wealthvar = "liquid_offshore_wealth_pareto",
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


Prop_offshoreweath <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "offshore_wealth",
        filters = filter_both, # age and non-negative
        descending = T #richest to poorest
      )
    })
  })
names(Prop_offshoreweath) <- allCountries

Prop_offshoreweath_liquid <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        wealthvar = "liquid_offshore_wealth",
        filters = filter_both, # age and non-negative
        descending = T #richest to poorest
      )
    })
  })
names(Prop_offshoreweath_liquid) <- allCountries

Prop_netwealth_asc <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt, 
        probs = seq(0,1,0.001),
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
        probs = seq(0,1,0.001),
        wealthvar = "liquid_assets",
        filters = filter_both, # age and non-negative
        descending = F #richest to poorest
      )
    })
  })
names(Prop_liquidassets_asc) <- allCountries


Prop_offshorewealth_asc <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        probs = seq(0,1,0.001),
        wealthvar = "offshore_wealth",
        filters = filter_both, # age and non-negative
        descending = F #richest to poorest
      )
    })
  })
names(Prop_offshorewealth_asc) <- allCountries

Prop_liquidoffshorewealth_asc <- allCountries %>%
  lapply(function(cnt) {
    mi_point(Wealthfiles, FUN = function(dset) {
      calc_WPercentile(
        dset,
        cntr=cnt,
        probs = seq(0,1,0.001),
        wealthvar = "liquid_offshore_wealth",
        filters = filter_both, # age and non-negative
        descending = F #richest to poorest
      )
    })
  })
names(Prop_liquidoffshorewealth_asc) <- allCountries

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

list(
  liquid_assets = Prop_liquidassets_asc$DE,
  net_wealth = Prop_netwealth_asc$DE,
  net_wealth_offshore = Prop_offshorewealth_asc$DE,
  liquid_assets_offshore = Prop_liquidoffshorewealth_asc$DE
) %>%
as_tibble %>%
lorenz_plotter(colr = c("gray","black"), titl = "Germany")

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


