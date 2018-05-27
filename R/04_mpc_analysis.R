source("R/00_libraries.R")
source("R/functions/dynprog.R")
source("R/functions/calibration.R")
source("R/classes/shocks_util_prod.R")
source("R/classes/krussell_smith.R")

library("reshape2")

# input values -----------------------------------------------------------

# Parameter values for psi and xi
# Obtained from (Carroll et al., 2014, Online Appendix, p. 6)
default_sigma_psi <- 0.01/4
default_sigma_xi <- 0.01*4
parameters <- list(
  AT = list(sigma_psi = default_sigma_psi, sigma_xi = default_sigma_xi,
            liq = c(beta_mid=0.9671509686, beta_range=0.04406199825),
            liq_off = c(beta_mid=0.953226385, beta_range=0.04243624967)),
  # BE = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  # CY = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  # DE = list(sigma_psi=default_sigma_psi, sigma_xi=0.05*4),
  # ES = list(sigma_psi=default_sigma_psi, sigma_xi=0.05*4),
  FI = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.953226385, beta_range=0.04243624967),
            liq_off = c(beta_mid=0.9497632451, beta_range=0.04361658526)),
  # FR = list(sigma_psi=default_sigma_psi, sigma_xi=0.031*4),
  GR = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.9445581602, beta_range=0.05528631474),
            liq_off = c(beta_mid=0.9041646851, beta_range=0.04576946483)),
  # IT = list(sigma_psi=default_sigma_psi, sigma_xi=0.075*4),
  # LU = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  # MT = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  # NL = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  PT = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi,
            liq = c(beta_mid=0.952499087, beta_range=0.04481387652),
            liq_off = c(beta_mid=0.9445581602, beta_range=0.05528631474))
  # SI = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  # SK = list(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi)
)
variables <- c("liq", "liq_off") # wealth variables to look at
models_file <- "data/generated/final_models.RData"
# the program will skip seaborn graphs if this variable is set to empty string ""
python_path <- "D:\\Run\\WinPython\\WinPython-64bit-3.6\\python-3.6.1.amd64\\python.exe"


# functions -------------------------------------------------------------

getFinalModels <- function(country_params) {
  models = list()
  for (wealth_var in variables) { # variables should be an argument
    estimated_betas = country_params[[wealth_var]]
    beta_min = estimated_betas["beta_mid"] - estimated_betas["beta_range"]
    beta_max = estimated_betas["beta_mid"] + estimated_betas["beta_range"]
    beta_n = 7 # we only consider the model with seven agent types
    
    # verify that estimated betas make sense
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


# solve on final parameters -----------------------------------------------

# only create models variable if it does not already exist and
# load from file if possible
if (!exists("models") & !file.exists(models_file)) {
  models <- lapply(parameters, getFinalModels)
  save(models, file = models_file)
} else if (!exists("models")) {
  load(models_file)
}

# caculate MPCs if they do not already exist
if (!exists("mpcs")) mpcs <- lapply(models, getMpcs)

# create one dataframe with all the MPCs
mpcs_wide <- melt(mpcs)
names(mpcs_wide) <- c("Estimate", "MPC", "Country")

# this object contains all the plots and summary statistics
mpc_analysis <- list()

# histograms
mpc_analysis$histograms <- lapply(mpcs, plotMpcHist)

# summary stats: mean, median, 
# TODO standard deviation
mpc_analysis$summaries <- lapply(mpcs, summary)

# box plots, v plots
mpc_analysis$boxplots <- lapply(mpcs, plotMpcBox)
mpc_analysis$vplots <- lapply(mpcs, plotMpcViolin)


# graphs with seaborn -----------------------------------------------------
if (python_path != "") {
  library(reticulate)
  use_python(python_path)
  sns <- import('seaborn')
  plt <- import('matplotlib.pyplot')
  pd <- import('pandas')
  
  # Seaborn set-up
  #sns$set(style="whitegrid")
  
  # strip plot
  # works but not a good idea - there are too many agents
  # stripplot <- sns$stripplot(x="MPC", y="Country", hue="Estimate", data=mpcs_wide, 
  #                            jitter=2, split=T, alpha=.25, zorder=1)
  # plt$show()
  
  # vplot
  mpc_analysis$vplot_py <- sns$violinplot(x="Country", y="MPC", hue="Estimate", data=mpcs_wide,
                                          split=TRUE)
  
  # plot the last graph
  #plt$show()
}

save(mpcs, mpcs_wide, mpc_analysis, file="data/generated/mpc_analysis.RData")