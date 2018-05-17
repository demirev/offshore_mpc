source("R/00_libraries.R")
source("R/functions/dynprog.R")
source("R/functions/calibration.R")
source("R/classes/shocks_util_prod.R")
source("R/classes/krussell_smith.R")

# input values -----------------------------------------------------------
Targets <- list(
  # Original CST data
  AT_cst = c(0, 78.8513, 91.7883, 96.6797, 99.0731, 100)/100,
  BE_cst = c(0, 62.7084, 81.8823, 92.2391, 97.9084, 100)/100,
  CY_cst = c(0, 72.098, 88.0831, 95.1412, 98.6535, 100)/100,
  DE_cst = c(0, 81.2646, 92.2557, 96.8536, 99.1951, 100)/100,
  ES_cst = c(0, 60.7636, 80.5163, 91.6185, 97.8128, 100)/100,
  FI_cst = c(0, 72.2088, 88.1598, 95.2109, 98.6785, 100)/100,
  FR_cst = c(0, 70.3265, 86.7258, 94.5155, 98.5333, 100)/100,
  GR_cst = c(0, 57.9778, 78.7277, 90.8523, 97.647, 100)/100,
  IT_cst = c(0, 63.4341, 82.1672, 92.3771, 98.0074, 100)/100,
  LU_cst = c(0, 70.2289, 86.8872, 94.6129, 98.5123, 100)/100,
  MT_cst = c(0, 63.005, 82.0733, 92.331, 97.9327, 100)/100,
  NL_cst = c(0, 67.8301, 85.3827, 93.9355, 98.3406, 100)/100,
  PT_cst = c(0, 65.3448, 83.6945, 93.1091, 98.1161, 100)/100,
  SI_cst = c(0, 57.9083, 78.7128, 90.8434, 97.6431, 100)/100,
  SK_cst = c(0, 54.9003, 76.8789, 90.1142, 97.5149, 100)/100
)

# Betas found by CST
betas_cst <- list(
  AT = c(beta_mid=0.988477, beta_range=0.0013),
  BE = c(beta_mid=0.98975, beta_range=0.0005),
  CY = c(beta_mid=0.989225, beta_range=0.0009),
  DE = c(beta_mid=0.986959, beta_range=0.0019),
  ES = c(beta_mid=0.9896, beta_range=0.0005),
  FI = c(beta_mid=0.9893, beta_range=0.0009),
  FR = c(beta_mid=0.98915, beta_range=0.0009),
  GR = c(beta_mid=0.9899, beta_range=0.0003),
  IT = c(beta_mid=0.989225, beta_range=0.0007),
  LU = c(beta_mid=0.989375, beta_range=0.0008),
  MT = c(beta_mid=0.98975, beta_range=0.0005),
  NL = c(beta_mid=0.9896, beta_range=0.0007),
  PT = c(beta_mid=0.98963, beta_range=0.0006),
  SI = c(beta_mid=0.9899, beta_range=0.0003),
  SK = c(beta_mid=0.989975, beta_range=-8.13152E-20)
)

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

#  Using the quantiles from CST ----
sink("calibration_checkpoints/log_IT20.txt")
calibrated_IT <- calibrate_genetic(
  FUN = function(betaPair) {
    run_calibration(beta_mid = betaPair[1], beta_range = betaPair[2], 
      sigma_xi = parameters$IT["sigma_psi"], 
      probs = seq(0, 1, 0.2) # original wealth data comes in quintiles
    )
  }, # wrapper around run_calibration that takes only 1 parameter
  lossF = lossKS(Targets$IT_cst), # function that will be used to evaluate
  individual_generator = generateKSParams(
    beta_mid_span = c(0.9, 0.99), 
    beta_rng_span = c(0.0001, 0.02)
  ), # function that will generate candidate parameters
  npop = 10, # size of population
  nsurvive = 4, # number of survivors per generation
  generations = 5, # number of generations to train
  tol = 1e-2, # will stop early if loss is less than this
  nparents = 2, # number of parents per children
  nchild = 4,
  checkpoint = "calibration_checkpoints/italy_20.csv", # file to write results to
  recordOutput = T # add quantiles to file
)
sink(NULL)

sink("calibration_checkpoints/log_ES20.txt")
calibrated_ES <- calibrate_genetic(
  FUN = function(betaPair) {
    run_calibration(beta_mid = betaPair[1], beta_range = betaPair[2], 
                    sigma_xi = parameters$ES["sigma_psi"], 
                    probs = seq(0, 1, 0.2) # original wealth data comes in quintiles
    )
  }, # wrapper around run_calibration that takes only 1 parameter
  lossF = lossKS(Targets$ES_cst), # function that will be used to evaluate
  individual_generator = generateKSParams(
    beta_mid_span = c(0.9, 0.99), 
    beta_rng_span = c(0.0001, 0.02)
  ), # function that will generate candidate parameters
  npop = 10, # size of population
  nsurvive = 4, # number of survivors per generation
  generations = 5, # number of generations to train
  tol = 1e-2, # will stop early if loss is less than this
  nparents = 2, # number of parents per children
  nchild = 4,
  checkpoint = "calibration_checkpoints/spain_20.csv", # file to write results to
  recordOutput = T # add quantiles to file
)
sink(NULL)

# solve on final parameters -----------------------------------------------


