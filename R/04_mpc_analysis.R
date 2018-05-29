source("R/00_libraries.R")
source("R/functions/dynprog.R")
source("R/functions/calibration.R")
source("R/classes/shocks_util_prod.R")
source("R/classes/krussell_smith.R")
source("R/functions/utils.R")
source("R/functions/descriptives.R")

library("reshape2")
library("ggrepel")

# input values -----------------------------------------------------------

# Parameter values 
# psi and xi are obtained from (Carroll et al., 2014, Online Appendix, p. 6)
default_sigma_psi <- 0.01/4
default_sigma_xi <- 0.01*4
# exclude CY, MT and LU which do not have offshore wealth estimates
parameters <- list(
  AT = list(sigma_psi = default_sigma_psi, sigma_xi = default_sigma_xi,
            liq = c(beta_mid=0.9615444847, beta_range=0.03845551433),
            liq_off = c(beta_mid=0.953226385, beta_range=0.04243624967)),
  BE = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.960384507, beta_range=0.03869966263),
            liq_off = c(beta_mid=0.9481308482, beta_range=0.04404311409)),
  # CY = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
  #           liq = c(beta_mid=0.9563485375, beta_range=0.04365146148),
  #           liq_off = c(beta_mid=0.9563485375, beta_range=0.04365146148)),
  DE = list(sigma_psi=default_sigma_psi, sigma_xi=0.05*4,
            liq = c(beta_mid=0.9692282535, beta_range=0.03077174549),
            liq_off = c(beta_mid=0.964282787, beta_range=0.03571721204)),
  ES = list(sigma_psi=default_sigma_psi, sigma_xi=0.05*4,
            liq = c(beta_mid=0.9480827236, beta_range=0.05191727536),  # x
            liq_off = c(beta_mid=0.9378592744, beta_range=0.06214072462)),
  FI = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.953226385, beta_range=0.04243624967),
            liq_off = c(beta_mid=0.9481308482, beta_range=0.04404311409)),
  FR = list(sigma_psi=default_sigma_psi, sigma_xi=0.031*4,
            liq = c(beta_mid=0.9473855263, beta_range=0.03758246963),
            liq_off = c(beta_mid=0.9453141982, beta_range=0.04461877977)),
  GR = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.9445581602, beta_range=0.05528631474),
            liq_off = c(beta_mid=0.9041646851, beta_range=0.04576946483)),
  IT = list(sigma_psi=default_sigma_psi, sigma_xi=0.075*4,
            liq = c(beta_mid=0.9328957257, beta_range=0.06404633582),
            liq_off = c(beta_mid=0.9322703618, beta_range=0.06687936708)),
  # LU = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
  #           liq = c(beta_mid=0.9607930198, beta_range=0.03147151368),
  #           liq_off = c(beta_mid=0.9607930198, beta_range=0.03147151368)),
  # MT = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
  #           liq = c(beta_mid=0.9445581602, beta_range=0.05528631474),
  #           liq_off = c(beta_mid=0.9445581602, beta_range=0.05528631474)),
  NL = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.9563485375, beta_range=0.04365146148),
            liq_off = c(beta_mid=0.9530021276, beta_range=0.04480030176)),
  PT = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.952499087, beta_range=0.04481387652),
            liq_off = c(beta_mid=0.9445581602, beta_range=0.05528631474)),
  SI = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.9445581602, beta_range=0.05528631474),
            liq_off = c(beta_mid=0.9445581602, beta_range=0.05528631474)),
  SK = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.952499087, beta_range=0.04481387652),
            liq_off = c(beta_mid=0.952499087, beta_range=0.04481387652))
)
variables <- c("liq", "liq_off") # wealth variables to look at
liq_label <- "Liquid Assets"
liq_off_label <- "Liquid Assets + Offshore Wealth"
models_file <- "data/generated/final_models.RData"
draw_gini_plot <- TRUE  # option to turn off gini plots because this takes for ever


# functions -------------------------------------------------------------

getFinalModels <- function(country_params) {
  models = list()
  for (wealth_var in variables) { # variables should be an argument
    estimated_betas = country_params[[wealth_var]]
    beta_min = estimated_betas["beta_mid"] - estimated_betas["beta_range"]
    beta_max = estimated_betas["beta_mid"] + estimated_betas["beta_range"]
    beta_n = 7 # we only consider the model with seven agent types
    
    # TODO verify that estimated betas make sense
    #if (beta_min <= 0) stop("lowest beta cannot be smaller than 0")
    #if (beta_max >= 1) stop("highest beta cannot be greater than 1")
    
    betas <- seq(
      from = max(beta_min, 0),
      to = min(beta_max, 1), 
      length.out = beta_n
    )
    
    models[[wealth_var]] <- KS_Economy$new(
      Agents = lapply(
        betas,
        function(beta) {
          AgentType$new(
            p = 1, # personal productivity initial
            l = 1/0.9, # fraction of labor supplied
            L = 1, # units of labor available
            psi = LogNormalShock$new(
              mu = 0, 
              sigma = country_params$sigma_psi
            ), # productivity permanent shock
            xi = EmploymentShock$new(
              mu = 0, 
              sigma = country_params$sigma_xi,
              mu_ins = 0.15, # percent insurance
              OM = 0.07 # probability of unemployment
            ), # productivity transitory shock
            utilf = Iso_Elastic$new(rho = 1), # utility function
            pol = function(m,k) return(0.5*m),
            D = 0.00625, # probability of death
            n = 3000, # number to simulate
            max_m = 35,
            num_out = 40,
            beta = beta # discount factor
          )
        }
      ), # list of agents
      P = 1, # aggregate productivity multiplier
      FF = Cobb_Douglas$new(alpha = 0.36), # production function
      PSI = NoShock$new(mu = 1), # aggregate permanent productivity schock
      XI  = NoShock$new(mu = 1), # aggregate transitory productivity schock
      delta = 0.025, # depreciation
      Z = 1 # TFP
    )
    
    # Solve the model
    models[[wealth_var]]$optimizeStohastic(
      k = 38.9, # initial value for capital
      tol = .02, # tolerance for the update of k_law
      tol_policy = .01, # tolerance for individual agent policies
      tol_ss = .02, # tolerance for steady state wealth distribution convergence
      lrate = .3, # speed of update for k_law
      max_m = 35, # discretization of m-space (highest gird point)
      num_out = 45, # discretization of m-space (number of grid points)
      fit_policy = fit_spline, # interpolation funciton
      verbose = 1, # print detailed messages or not
      probs = seq(0,1,0.1) # output percentiles - make sure they match Target
    )
  }
  models
}

getMpcs <- function(models) {
  "Returns a dataframe with the mpc of every individual in each calibrated model"
  mpcs <- lapply(models, function(model){
    unlist(model$calcMPC()[["mpc_list"]])
  })
  mpcs <- as.data.frame(mpcs)
}

widenMpcList <- function(mpcs) {
  "Puts mpcs into a wide data frame"
  mpcs <- melt(mpcs, variable.name="Estimate")
  names(mpcs)[names(mpcs) == 'value'] <- 'MPC'
  mpcs
}

plotMpcHist <- function(mpcs_list, country_name="NA", alpha=0.25, binwidth=.1) {
  "Plots a joint histogram of mpcs for multiple models"
  
  mpcs_list <- widenMpcList(mpcs_list)
  
  ggplot(mpcs_list, aes(x=MPC, fill = Estimate)) + 
    geom_histogram(alpha=alpha, binwidth=binwidth, position = "identity") # overlapping histogram
}

plotMpcViolin <- function(mpcs_list, country_name="NA") {
  "Plots a joint violin plot of mpcs for multiple models"
  
  mpcs_list <- widenMpcList(mpcs_list)
  
  ggplot(mpcs_list, aes(y=MPC, x = Estimate)) + 
    geom_violin()
}

plotMpcBox <- function(mpcs_list, country_name="NA") {
  "Plots a joint box plot of mpcs for multiple models"
  
  mpcs_list <- widenMpcList(mpcs_list)
  
  ggplot(mpcs_list, aes(y=MPC, x = Estimate)) + 
    geom_boxplot()
}

plotMpcVsGini <- function(mpcs_wide, gini_liq, gini_liq_off, liq_label, liq_off_label) {
  "Plots average MPCs against Ginis for all countries"
  df <- aggregate(mpcs_wide[,2], mpcs_wide[,c(1,3)], FUN = mean)  # mean mpcs
  names(df)[[3]] <- "AvgMPC"
  
  # combine two gini arrays & match format of mpc_wide
  gini_df <- as.data.frame(cbind(gini_liq, gini_liq_off))
  names(gini_df) <- c("liq", "liq_off")  # match Estimate values in mpcs
  gini_df$Country <- rownames(gini_df)  # create country column
  gini_df <- melt(gini_df, id.vars = "Country")  # match columns of mpcs
  names(gini_df) <- c("Country", "Estimate", "GINI")  # match column names
  gini_df$Estimate <- as.factor(gini_df$Estimate)  # match Estimate format
  levels(gini_df$Estimate)[levels(gini_df$Estimate)=="liq"] <- liq_label
  levels(gini_df$Estimate)[levels(gini_df$Estimate)=="liq_off"] <- liq_off_label
  
  # merge
  df <- merge(df, gini_df, by=c("Country", "Estimate"), all = TRUE)
  
  # plot the whole thing
  ggplot(df, aes(x=GINI, y=AvgMPC, color=Estimate, shape=Estimate)) + 
    geom_point() + 
    geom_label_repel(  # instead of geom_text(aes(label = Country))
      aes(label = Country),
      size = 3,
      show.legend = FALSE,
      label.size = 0  # supposedly deactivates the border of the labels
    ) +
    theme_bw()  # remove background
}


# solve and analyse models ----------------------------------------------------

# only create models variable if it does not already exist and
# load from file if possible
if (!exists("models") & !file.exists(models_file)) {
  models <- lapply(parameters, getFinalModels)
  save(models, file = models_file)
} else if (!exists("models")) {
  load(models_file)
}

# caculate MPCs if they do not already exist
mpcs <- lapply(models, getMpcs)

# create one dataframe with all the MPCs
mpcs_wide <- melt(mpcs)
names(mpcs_wide) <- c("Estimate", "MPC", "Country")
# give the Estimate factor readable values
levels(mpcs_wide$Estimate)[levels(mpcs_wide$Estimate)=="liq"] <- liq_label
levels(mpcs_wide$Estimate)[levels(mpcs_wide$Estimate)=="liq_off"] <- liq_off_label

# this object contains all the plots and summary statistics
mpc_analysis <- list()

# histograms
#mpc_analysis$histograms <- lapply(mpcs, plotMpcHist)

# summary stats: mean, median, 
# TODO standard deviation?
mpc_analysis$summaries <- lapply(mpcs, summary)

# box plots, v plots
#mpc_analysis$boxplots <- lapply(mpcs, plotMpcBox)
#mpc_analysis$vplots <- lapply(mpcs, plotMpcViolin)

# MPCs versus GINI
if (draw_gini_plot) {
  load("data/generated/wealthfiles_pareto_Christoph.RData")
  if (!exists("Wealthfiles")) Wealthfiles <- bigImport()
  allCountries <- unique(Wealthfiles[[1]]$country)
  if (!exists("gini_liq")) gini_liq <- allCountries %>%
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
  if (!exists("gini_liq_off")) gini_liq_off <- allCountries %>%
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
  mpc_analysis$mpc_vs_gini <- plotMpcVsGini(mpcs_wide, gini_liq, gini_liq_off,
                                            liq_label, liq_off_label)
}


# wealth distributions by Joel ----------------------------------------------------

load('data/generated/decile_targets.RData')
#We want to loop only through the countries we have already calibrated
countries <- as.list(names(models))
names(countries) <- countries
#Create list containing target and calibrated distributions and their difference
distributions <- list()
for (i in countries) {
  #Non-offshore distributions
  distributions[[i]]$liq$Target <- Targets$liq[i,]
  distributions[[i]]$liq$Calibration <- models[[i]]$liq$quantileSummary()$quantiles
  distributions[[i]]$liq$Difference <- distributions[[i]]$liq$Target -
    distributions[[i]]$liq$Calibration
  #Offshore distributions
  distributions[[i]]$liq_offshore$Target <- Targets$liq_offshore[i,]
  distributions[[i]]$liq_offshore$Calibration <- models[[i]]$liq_off$quantileSummary()$quantiles
  distributions[[i]]$liq_offshore$Difference <- distributions[[i]]$liq_offshore$Target -
    distributions[[i]]$liq_offshore$Calibration
}
names_distr <- rep('', length(countries)*6)
seq <- seq_along(countries)*5-5
for (i in 1:length(countries)) {
  names_distr[i+seq[i]] <- countries[i]
}
Adjustment <- rep(c('Without Offshore', '','','With Offshore','',''),length(countries))
distr_type <- rep(c('Target','Calibrated','Difference'), length(countries))
df <- t(data.frame(distributions))
df1 <- data.frame(Adjustment, distr_type, df)
df2 <- as.matrix(df1)
rownames(df2) <- names_distr 
colnames(df2) <- c('Adjusmtent', 'Distribution', colnames(df))
distributions_mat <-apply(df2, c(1,2), paste, collapse = "")
rm(df, df1, df2)


save(mpcs, mpcs_wide, mpc_analysis, distributions, distributions_mat,
     file="data/generated/mpc_analysis.RData")
