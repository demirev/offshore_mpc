source("R/00_libraries.R")
source("R/functions/dynprog.R")
source("R/functions/calibration.R")
source("R/classes/shocks_util_prod.R")
source("R/classes/krussell_smith.R")

# input target values -----------------------------------------------------
Targets <- list(
  AT = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
)

# define parameters -------------------------------------------------------
# the below can be of-course written much more concisely, but for purposes
# of keeping all used parameters transperent everything is spelled out

run_calibration <- function(beta_mid, beta_range, beta_n = 7) {
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
            sigma = 0.025
          ), # productivity permanent shock
          xi = EmploymentShock$new(
            mu = 0, 
            sigma = 0.04,
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
    probs = seq(0,1,0.1) # output percentiles - make sure they match Target
  )
  
  return(optResult)
}

# run ---------------------------------------------------------------------
calibrated_AT <- calibrate_genetic(
  FUN = function(betaPair) {
    run_calibration(beta_mid = betaPair[1], beta_range = betaPair[2])
  }, # wrapper around run_calibration that takes only 1 parameter
  lossF = lossKS(Targets$AT), # function that will be used to evaluate
  individual_generator = generateKSParams(
    beta_mid_span = c(0.9, 0.99), 
    beta_rng_span = c(0.01, 0.05)
  ), # function that will generate candidate parameters
  npop = 4, # size of population
  nsurvive = 2, # number of survivors per generation
  generations = 3, # number of generations to train
  tol = 1e-2, # will stop early if loss is less than this
  nparents = 2, # number of parents per children
  nchild = 1, # number of children to spawn per generation
  initial_pop = NULL, # starting population can be supplied directly
  checkpoint = "calibration_checkpoints/test.csv", # file to write results to
  recordOutput = T # add quantiles to file
)


# solve on final parameters -----------------------------------------------


