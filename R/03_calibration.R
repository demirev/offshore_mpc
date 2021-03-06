source("R/00_libraries.R")
source("R/functions/dynprog.R")
source("R/functions/calibration.R")
source("R/classes/shocks_util_prod.R")
source("R/classes/krussell_smith.R")

# input values -----------------------------------------------------------

# Country
country <- "NL"

# Load Targets from a file, since it is a rather long list
load('data/generated/decile_targets.RData')

# Parameter values
# Obtained from (Carroll et al., 2014, Online Appendix, p. 6)
default_sigma_psi = 0.01/4
default_sigma_xi = 0.01*4
parameters <- list(
  AT = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  BE = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  CY = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  DE = c(sigma_psi=default_sigma_psi, sigma_xi=0.05*4),
  ES = c(sigma_psi=default_sigma_psi, sigma_xi=0.05*4),
  FI = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  FR = c(sigma_psi=default_sigma_psi, sigma_xi=0.031*4),
  GR = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  IT = c(sigma_psi=default_sigma_psi, sigma_xi=0.075*4),
  LU = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  MT = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  NL = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  PT = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  SI = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi),
  SK = c(sigma_psi=default_sigma_psi, sigma_xi=default_sigma_xi)
)


# define parameters -------------------------------------------------------
run_calibration <- function(beta_mid, beta_range, beta_n = 7, 
  sigma_psi = 0.01/4, # variance of log permanent shocks
  sigma_xi = 0.01*4, # variance of log transitory shocks
  probs = seq(0,1,0.1) # output percentiles - make sure they match Target
) {
  betas <- seq(
    from = max(beta_mid - beta_range, 0),
    to = min(beta_mid + beta_range, 1), 
    length.out = beta_n
  )
  
  model <- KS_Economy$new(
    Agents = lapply(
      betas,
      function(beta) {
        AgentType$new(
          p = 1, # personal productivity initial
          l = 1/0.9, # fraction of labor supplied
          L = 1, # units of labor available
          psi = LogNormalShock$new(
            mu = 0, 
            sigma = sigma_psi
          ), # productivity permanent shock
          xi = EmploymentShock$new(
            mu = 0, 
            sigma = sigma_xi,
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
  
  optResult <- model$optimizeStohastic(
    k = 38.9, # initial value for capital
    tol = .02, # tolerance for the update of k_law
    tol_policy = .01, # tolerance for individual agent policies
    tol_ss = .02, # tolerance for steady state wealth distribution convergence
    lrate = .3, # speed of update for k_law
    max_m = 35, # discretization of m-space (highest gird point)
    num_out = 45, # discretization of m-space (number of grid points)
    fit_policy = fit_spline, # interpolation funciton
    verbose = 1, # print detailed messages or not
    probs = probs # output percentiles - make sure they match Target
  )
  
  return(optResult)
}

# run ---------------------------------------------------------------------

# Liquid assets, without offshore ----
date = format(Sys.time(), "%y%m%d")
calibrated_liq <- calibrate_genetic(
  FUN = function(betaPair) {
    run_calibration(
      beta_mid = betaPair[1], beta_range = betaPair[2], 
      sigma_xi = get(country, parameters)["sigma_xi"], 
      probs = seq(0, 1, 0.1) # deciles 
    )
  }, # wrapper around run_calibration that takes only 1 parameter
  lossF = lossKS(Targets$liq[country,]), # function that will be used to evaluate
  individual_generator = generateKSParams(
    beta_mid_span = c(0.90, 0.98), 
    beta_rng_span = c(0.01, 0.05)
  ), # function that will generate candidate parameters
  npop = 10, # size of population
  nsurvive = 4, # number of survivors per generation
  generations = 5, # number of generations to train
  tol = 5e-2, # will stop early if loss is less than this
  nparents = 3, # number of parents per children
  nchild = 2,
  checkpoint = str_glue("calibration_checkpoints/{country}_liq.csv"), # file to write results to
  recordOutput = T, # add quantiles to file
  logMessages = T,
  log_file = str_glue("calibration_checkpoints/log_{country}_liq_{date}.txt")
)

# Liquid assets, with offshore ----
date = format(Sys.time(), "%y%m%d")
calibrated_liq_off <- calibrate_genetic(
  FUN = function(betaPair) {
    run_calibration(
      beta_mid = betaPair[1], beta_range = betaPair[2],
      sigma_xi = get(country, parameters)["sigma_xi"],
      probs = seq(0, 1, 0.1) # deciles
    )
  }, # wrapper around run_calibration that takes only 1 parameter
  lossF = lossKS(Targets$liq_offshore[country,]), # function that will be used to evaluate
  individual_generator = generateKSParams(
    beta_mid_span = c(0.90, 0.98),
    beta_rng_span = c(0.01, 0.05)
  ), # function that will generate candidate parameters
  npop = 10, # size of population
  nsurvive = 4, # number of survivors per generation
  generations = 5, # number of generations to train
  tol = 5e-2, # will stop early if loss is less than this
  nparents = 3, # number of parents per children
  nchild = 2,
  checkpoint = str_glue("calibration_checkpoints/{country}_liq_off.csv"), # file to write results to
  recordOutput = T, # add quantiles to file
  logMessages = T,
  log_file = str_glue("calibration_checkpoints/log_{country}_liq_off_{date}.txt")
)

