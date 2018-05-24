source("R/00_libraries.R")
source("R/functions/dynprog.R")
source("R/functions/calibration.R")
source("R/classes/shocks_util_prod.R")
source("R/classes/krussell_smith.R")

# input values -----------------------------------------------------------

# Parameter values for psi and xi
# Obtained from (Carroll et al., 2014, Online Appendix, p. 6)
default_sigma_psi = 0.01/4
default_sigma_xi = 0.01*4
parameters <- list(
  AT = list(sigma_psi = default_sigma_psi, sigma_xi = default_sigma_xi,
         liq = c(beta_mid=0.9671509686, beta_range=0.04406199825),
         liq_off = c(beta_mid=0.953226385, beta_range=0.04243624967))
)


# solve on final parameters -----------------------------------------------

# contains one model for each country and model type, e.g. models$AT$betas_liq
models = list() 

country = "AT"
probs = seq(0,1,0.1)

# everythin below this lines should go in a country loop 
models[[country]] = list()

for (specification in c("liq", "liq_off")) {
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
  
  models[[country]][[specification]] <- KS_Economy$new(
    Agents = lapply(
      betas,
      function(beta) {
        AgentType$new(
          p = 1, # personal productivity initial
          l = 1/0.9, # fraction of labor supplied
          L = 1, # units of labor available
          psi = LogNormalShock$new(
            mu = 0, 
            sigma = get(country, parameters)$sigma_psi
          ), # productivity permanent shock
          xi = EmploymentShock$new(
            mu = 0, 
            sigma = get(country, parameters)$sigma_xi,
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
  models[[country]][[specification]]$optimizeStohastic(
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
}

save(models, file = "models.RData")
