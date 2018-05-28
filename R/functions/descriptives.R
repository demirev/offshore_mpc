library(tidyr)
library(DescTools) # for Gini
library(plotly)
library(kableExtra)
library(RColorBrewer)

# Import ------------------------------------------------------------------
bigImport <- function(targetdir = "data/HFCS_UDB_1_3_ASCII/") {
  
  # just wraps a bunch of data prep code - will be passed to rmd files
  
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
    bonds = "DA2103",
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
  
  # import off-shore wealth figures
  # temp - should be read from original Zucman file in final version
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
    lapply(function(x) x*0.773477635) %>% # to eur 2009
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
            shares + managed_accounts + life_insurance + bonds,
          weight = household_weight,
          permanent_income_flag = permanent_income_flag
        ) %>%
        # mutate(
        #   liquid_assets = liquid_assets * (liquid_assets > 0),
        #   net_wealth = net_wealth * (net_wealth > 0)
        # ) %>%
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
          ),
          offshore_wealth_pareto = calc_offshore(
            wealthvar = net_wealth,#liquid_assets, 
            offshore = OffShoreWealth[country][[1]],
            weight = weight/10
          ),
          liquid_offshore_wealth_pareto = calc_offshore(
            wealthvar = liquid_assets, 
            offshore = OffShoreWealth[country][[1]],
            weight = weight/10
          )
        ) %>%
        ungroup
    })
}


# Off-shore Wealth Merge --------------------------------------------------
calc_offshore <- function(wealthvar, offshore, 
                          probs = seq(0,1,length.out = length(wealthvar)), 
                          weight,
                          share_bottom_90 = 1.58226, coeff = 2.7672) {
  
  weight[weight < 1] <- 1 # introduces a bug if weight = 0 after rounding
  
  # create weighted vector
  wealthvar_weighted <- rep(wealthvar, round(weight))
  original_id <- rep(1:length(wealthvar), round(weight))
  print(pryr::object_size(wealthvar_weighted))
  
  # calculate quantiles of weighted vector
  qs <- quantile(wealthvar, seq(0,1, length.out = length(wealthvar_weighted)))
  ps <- qs %>% names %>% str_replace("%","") %>% as.numeric
  
  # keep track of ordering so we can rearrange
  original_order <- original_id[order(wealthvar_weighted)]
  wealthvar_ordered <- wealthvar_weighted[order(wealthvar_weighted)]
  
  
  cumshares = vector("numeric",length(wealthvar_ordered))
  
  # split equally among first 90 percent
  cumshares[wealthvar_ordered <= max(qs[ps <= 90])] <- cumsum(
    rep(
      share_bottom_90/sum(wealthvar_ordered <= max(qs[ps <= 90])), 
      sum(wealthvar_ordered <= max(qs[ps <= 90]))
    )
  )
  
  # top 10 percent as per coefficients of a pareto model
  cumshares[wealthvar_ordered > max(qs[ps <= 90])] <- share_bottom_90 /
    (1 - ps[wealthvar_ordered > max(qs[ps <= 90])]/100)^(1/coeff) 
  
  # handle issue with multiple individuals at 100% - approach 1
  # cumshares[cumshares > 100] <- 100
  # cumshares <- cumshares/100
  # shares <- cumshares - lag(cumshares)
  # shares[1] <- 0
  # shares[shares < 0] <- 0
  # shares[cumshares == 1] <- shares[shares != 0][length(shares[shares != 0])]/sum(cumshares == 1)
  
  # approach 2 - rescaling the tail
  above_90 <- wealthvar_ordered > max(qs[ps <= 90])
  cumshares[above_90] <- cumshares[above_90] - share_bottom_90
  cumshare_adj <- (cumshares[length(cumshares)]-share_bottom_90)/(100-share_bottom_90)
  cumshares[above_90] <- cumshares[above_90]/cumshare_adj
  cumshares[above_90] <- cumshares[above_90] + share_bottom_90
  cumshares <- cumshares/100
  shares <- cumshares - lag(cumshares)
  
  offshore_wealth <- wealthvar_ordered + shares*offshore
  
  # reorganize in original state
  results <- tibble(offshore_wealth = offshore_wealth, original_order = original_order) %>%
    group_by(original_order) %>%
    summarise(offshore = mean(offshore_wealth))
  
  return(results$offshore)
}

calc_offshore2 <- function(
  wealthvar, offshore, 
  weight,
  share_bottom_90 = 0.0158226,
  share_top_10 = 0.009891366,
  share_top_5 = 0.038324026,
  share_top_1 = 0.034766489,
  share_top_05 = 0.132246103,
  share_top_01 = 0.253311256,
  share_top_001 = 0.515638502
) {
  #browser()
  probs = c(0.9,0.95,0.99,0.995,0.999,0.9999)
  qs <- quantile(rep(wealthvar,weight), probs)
  
  # quick version - manual input of a lot of things
  cond1 <- wealthvar < qs["90%"]
  cond2 <- wealthvar >= qs["90%"] & wealthvar < qs["95%"]
  cond3 <- wealthvar >= qs["95%"] & wealthvar < qs["99%"]
  cond4 <- wealthvar >= qs["99%"] & wealthvar < qs["99.5%"]
  cond5 <- wealthvar >= qs["99.5%"] & wealthvar < qs["99.9%"]
  cond6 <- wealthvar >= qs["99.9%"] & wealthvar < qs["99.99%"]
  cond7 <- wealthvar >= qs["99.99%"]
  
  offshore_dist <- cond1*share_bottom_90 + cond2*share_top_10 + 
    cond3*share_top_5 + cond4*share_top_1 + cond5*share_top_05 +
    cond6*share_top_01 + cond7*share_top_001
  bygroupweight <- cond1*sum(weight[cond1]) + cond2*sum(weight[cond2]) +
    cond3*sum(weight[cond3]) + cond4*sum(weight[cond4]) +
    cond5*sum(weight[cond5]) + cond6*sum(weight[cond6]) +
    cond7*sum(weight[cond7])
  
  result <- wealthvar + offshore*offshore_dist/bygroupweight
  result
}

# Filters -----------------------------------------------------------------
filter_age <- function(dset, wealthvar) {
  dset %>%
    filter((age >= 25 & age <= 60) |
             (country == "MT" & age_bin >= 6 & age_bin < 13))
}

filter_income <- function(dset, wealthvar) {
  dset %>%
    filter(permanent_income_flag == 2 | country %in% c("FR","FI")) 
}

filter_negative <- function(dset, wealthvar) {
  dset[dset[[wealthvar]] >= 0, ]
}

filter_both <- function(dset, wealthvar) {
  dset %>%
    filter_age(wealthvar) %>%
    filter_negative(wealthvar)
}

filter_three <- function(dset, wealthvar) {
  dset %>%
    filter_both(wealthvar) %>%
    filter_income(wealthvar)
}

omni_filt <- function(dset, cntr = NULL, 
                      filters = function(dset, ...){return(dset)}, 
                      descending = T,
                      wealthvar = "net_wealth",
                      weightvar = "weight", cntryvar = "country") {
  # some boilerplate that goes in several places below
  dset <- filters(dset, wealthvar)
  if (!is.null(cntr)) {
    dset <- dset[dset[cntryvar] == cntr, ]
  } 
  dset <- dset[order(dset[[wealthvar]], decreasing = descending), ]
  return(dset)
}


# Distributions -----------------------------------------------------------
flag_quantile <- function(var, quantiles) {
  flags <- lapply(quantiles, function(q) var <= q) %>%
    as.data.frame %>%
    rowSums 
  
  # this works for one clump only :(
  notin <- (length(quantiles):1)[!(length(quantiles):1 %in% flags)]
  
  if (length(notin)) {
   
    isin  <- max(notin)+1
    dividers <- c(isin, notin)
    if (!flags[1] == length(quantiles)) dividers <- rev(dividers)
    oldlev <- flags == isin
    
    todivide <- length(flags[flags == isin])
    shares <- (todivide/(length(notin)+1)) %>% round
    for (newlev in dividers) {
      frm <- (which(c(isin, notin)==newlev)-1)*shares+1
      to  <- min(which(c(isin, notin)==newlev)*shares,todivide)
      flags[oldlev][frm:to] <- newlev
    }
     
  }
  
  flags <- as.factor(flags)
  levels(flags) <- rev(names(quantiles))
  
  flags
}

calc_Totals <- function(
  dset, 
  cntr = NULL,
  filters = function(dset,...){return(dset)}, 
  descending = T,
  wealthvar = "net_wealth", # for filtering
  sumvars = c("net_wealth","liquid_assets","offshore_wealth"),
  weightvar = "weight", cntryvar = "country"
) {
  # weighted column sums
  dset <- omni_filt(dset, cntr, filters, 
                    descending, wealthvar, weightvar, cntryvar)
  
  sapply(sumvars, function(sumvar) sum(dset[[sumvar]]*dset[[weightvar]]))
}

calc_Gini <- function(
  dset, 
  cntr = NULL,
  descending = T,
  filters = function(dset,...){return(dset)}, 
  wealthvar = "net_wealth",
  weightvar = "weight", cntryvar = "country"
) {
  # wrapper aroung Gini that does some filtering
  dset <- omni_filt(dset, cntr, filters, 
                    descending, wealthvar, weightvar, cntryvar)
  
  wealth <- dset[[wealthvar]]
  weight <- dset[[weightvar]]
  
  Gini(wealth, n = round(weight)) # point estimate
}


calc_WPercentile <- function(
  dset, cntr = NULL,  descending = T,
  probs = seq(0,1,0.05), cumulative = T,
  filters = function(dset,...){return(dset)},
  wealthvar = "net_wealth",
  weightvar = "weight", 
  cntryvar = "country",
  negative_to_zero = F
) {
  # weighted cumulative wealth distribution quantiles
  dset <- omni_filt(dset, cntr, filters, 
                    descending, wealthvar, weightvar, cntryvar)
  
  wealth <- rep(dset[[wealthvar]], round(dset[[weightvar]])) # weighted
  
  if (negative_to_zero){
    wealth[wealth < 0] = 0
  }
  
  if (cumulative) {
    quantile(x = cumsum(wealth), probs = probs) / sum(wealth)
  } else {
    quantile(x = wealth, probs = probs)
  }
  
}

calc_WtoI <- function(
  dset, cntr = NULL, descending = F,
  probs = seq(0,1,0.05), cumulative = T,
  filters = function(dset,...){return(dset)},
  high_cutoff = 10000, # value in C paper
  wealthvar = "net_wealth",
  incomevar = "income",
  weightvar = "weight",
  cntryvar = "country"
) {
  
  # wealth to income ratios
  # EG: Poorest 10% have 2 w/i ratio, poorest 20% have ... etc
  if (!is.null(high_cutoff)) {
    filt_high <- function(dset) {
      toohigh <- (dset$qwealth/dset$qincome > high_cutoff)
      dset$qwealth[toohigh] <- 10000 # cap
      dset$qincome[toohigh] <- 1
      dset
    }
  } else {
    filt_high <- function(dset) dset
  }
  
  dset <- omni_filt(dset, cntr, filters, 
                    descending, wealthvar, weightvar, cntryvar)
  
  wealth <- rep(dset[[wealthvar]], round(dset[[weightvar]])) # weighted
  income <- rep(dset[[incomevar]], round(dset[[weightvar]]))
  
  qs <- flag_quantile(wealth, quantile(x = wealth, probs = probs))
  
  if (cumulative) { # cumullative per quantile
    r <- tibble(
      ratio = cumsum(wealth)/(cumsum(income)/4), 
      qfl = qs
    ) %>%
      group_by(qfl) %>%
      filter(row_number()==ifelse(descending, n(),1))
  } else { # average per quantile
    r <- tibble(
      qwealth = wealth,
      qincome = (income/4),
      qfl = qs
    ) %>%
      filt_high() %>% # filter high wealth/income rations
      group_by(qfl) %>%
      summarise(ratio = sum(qwealth)/sum(qincome))
  }
  
  result <- r$ratio
  names(result) <- r$qfl
  return(result)
}


# Visulizations -----------------------------------------------------------
lorenz_plotter <- function(tib, titl = "Lorenz Curves", 
                           colr = c("gray","black")) {
  
  pallt <- colorRampPalette(col=colr)(ncol(tib))
  
  p <- tib %>%
    plot_ly(
      x = seq(0,1,length.out = nrow(tib)), 
      y = seq(0,1,length.out = nrow(tib)), name = 'baseline', 
      type = 'scatter', mode = 'lines',
      line = list(color = 'rgb(205, 12, 24)', width = 1)
    )
  for (col in seq_along(tib)) {
    p <- add_trace(
      p, y = tib[[col]], name = colnames(tib)[col],
      line = list(
        color = pallt[col],
        width = 2
      )             
    )
  }
  p <- layout(
    p, 
    title = titl,
    xaxis = list(
      title = "Cumulative % of Population",
      showgrid = F
    ),
    yaxis = list(
      title = "Cumulative % of Wealth",
      showgrid = F
    )
  )
  return(p)
}

bar_plotter <- function(tib, lbl_names, titl = "Gini Coefficients",
                        colr = c("seashell","seashell3")) {
  
  pallt <- colorRampPalette(col=colr)(ncol(tib))
  
  p <- tib %>%
    plot_ly(x = ~lbl_names)
  
  for (col in seq_along(tib)) {
    p <- add_trace(
      p, y = tib[[col]], type = 'bar',
      name = colnames(tib)[[col]],
      text = round(tib[[col]]*100, 2), textposition = 'auto',
      marker = list(
        color = pallt[col],
        line = list(
          #color = 'rgb(195,1,1)', 
          width = 1.5
        )
      )
    )
  }
  
  p <- layout(
    p,
    title = titl,
    barmode = 'group',
    xaxis = list(
      title = ""
    ),
    yaxis = list(
      title = "",
      showgrid = F
    )
  )
  
  return(p)
}

box_plotter <- function(tib, titl = "Wealth/Income Distribution", 
                        colr = c("gray","black")) {
  
  pallt <- colorRampPalette(col=colr)(length(tib))
  
  gen_box <- function(mn, p25, p50, p75, mx) {
    mn <- max(mn,p25-(p75-p25)*1.5)
    mx <- min(mx,p75+(p75-p25)*1.5)
    c(rep(mn,24), rep(p25,24),rep(p50,4),rep(p75,24),rep(mx,24))
  }
  
  gen_box_from_table <- function(tablevar) {
    nms <- c("0%","25%","50%","75%","100%")
    if (!all(nms %in% names(tablevar))) stop("Missing Percentiles")
    gen_box(tablevar[nms[1]],tablevar[nms[2]],tablevar[nms[2]],
            tablevar[nms[4]],tablevar[nms[5]])
  }
  
  p <- tibble(
    ratio = Reduce(c,lapply(tib, gen_box_from_table)),
    cntry = rep(names(tib), each = 100)
  ) %>%
    plot_ly(y = ~ratio, color = ~cntry, type = "box", colors = pallt) %>%
    layout(
      title = titl,
      barmode = 'group',
      autosize = F, width = 1000,
      xaxis = list(
        title = "Country", range = c(0,95)
      ),
      yaxis = list(
        title = "Wealth-to-Income Ratio",
        showgrid = F
      ),
      legend = list(
        orientation = 'h'
      )
    )
  
  return(p)
}
