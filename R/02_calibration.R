source("R/functions/dynprog.R")
source("R/functions/dynprog.R")
source("R/classes/shocks_util_prod.R")
source("R/classes/krussell_smith.R")


# input target values -----------------------------------------------------


# define parameters -------------------------------------------------------
# the below can be of-course written much more concisely, but for purposes
# of keeping all used parameters transperent everything is spelled out

run_calibration <- function(beta_mid, beta_range, beta_n) {
  betas <- seq(
    from = beta_mid - range, 
    to = beta_mid + range, 
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
          n = 10000, # number to simulate
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
}

# run ---------------------------------------------------------------------



# solve on final parameters -----------------------------------------------


