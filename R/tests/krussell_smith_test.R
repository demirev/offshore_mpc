source("R/00_libraries.R")
source("R/classes/shocks_util_prod.R")
source("R/classes/krussell_smith.R")
source("R/functions/dynprog.R")

kStart = 52.95991 # 38.9
nAgent = 5000

#----------------------------------------------------------------------
testKs <- KS_Economy$new(
  Agents = list(
    AgentType$new(
      M = kStart,
      p = 1,
      l = 1/0.9,
      L = 1,
      psi = LogNormalShock$new(sigma = 0.025),
      xi = EmploymentShock$new(sigma = 0.04),
      utilf = Iso_Elastic$new(rho = 1),
      pol = function(m, k) return(m*0.5),
      QTable = NULL,
      D = 0.00625,
      beta = 0.985,
      n = nAgent
    )
  ),
  K = kStart * (nAgent/0.9),
  P = 1,
  FF = Cobb_Douglas$new(alpha = 0.36),
  PSI = NoShock$new(mu = 1),
  XI  = NoShock$new(mu = 1),
  delta = 0.025
)
# 
testKs$stohastic_optimization(k = kStart) #52.95991
# testKs$update_policies(k = testKs$K / testKs$L)
# testKs$wealth_ss(verbose = T)

#----------------------------------------------------------------------

testKs2 <- KS_Economy$new(
  Agents = list(
    AgentType$new(
      M = kStart,
      p = 1,
      l = 1/0.9,
      L = 1,
      psi = LogNormalShock$new(sigma = 0.025),
      xi = EmploymentShock$new(sigma = 0.04),
      utilf = Iso_Elastic$new(rho = 1),
      pol = function(m, k) return(m*0.5),
      QTable = NULL,
      D = 0.00625,
      beta = 0.985,
      n = nAgent
    ),
    AgentType$new(
      M = kStart,
      p = 1,
      l = 1/0.9,
      L = 1,
      psi = LogNormalShock$new(sigma = 0.025),
      xi = EmploymentShock$new(sigma = 0.04),
      utilf = Iso_Elastic$new(rho = 1),
      pol = function(m, k) return(m*0.5),
      QTable = NULL,
      D = 0.00625,
      beta = 0.965,
      n = nAgent
    )
  ),
  K = kStart * (nAgent/0.9) * 2,
  P = 1,
  FF = Cobb_Douglas$new(alpha = 0.36),
  PSI = NoShock$new(mu = 1),
  XI  = NoShock$new(mu = 1),
  delta = 0.025
)
testKs2$updatePolicies(k = testKs2$K/testKs2$L, law_k = function(k) k, tol = 0.04)
testKs2$plotDist()
testKs2$Agents[[1]]$plotPolicy()
testKs2$Agents[[2]]$plotPolicy()
testKs2$findWealthSS(verbose = T)

#----------------------------------------------------------------------

testKs3 <- KS_Economy$new(
  Agents = list(
    AgentType$new(
      M = kStart,
      p = 1,
      l = 1/0.9,
      L = 1,
      psi = LogNormalShock$new(sigma = 0.025),
      xi = EmploymentShock$new(sigma = 0.04),
      utilf = Iso_Elastic$new(rho = 1),
      pol = function(m, k) return(m*0.5),
      QTable = NULL,
      D = 0.00625,
      beta = 0.985,
      n = nAgent
    ),
    AgentType$new(
      M = kStart,
      p = 1,
      l = 1/0.9,
      L = 1,
      psi = LogNormalShock$new(sigma = 0.025),
      xi = EmploymentShock$new(sigma = 0.04),
      utilf = Iso_Elastic$new(rho = 1),
      pol = function(m, k) return(m*0.5),
      QTable = NULL,
      D = 0.00625,
      beta = 0.975,
      n = nAgent
    ),
    AgentType$new(
      M = kStart,
      p = 1,
      l = 1/0.9,
      L = 1,
      psi = LogNormalShock$new(sigma = 0.025),
      xi = EmploymentShock$new(sigma = 0.04),
      utilf = Iso_Elastic$new(rho = 1),
      pol = function(m, k) return(m*0.5),
      QTable = NULL,
      D = 0.00625,
      beta = 0.980,
      n = nAgent
    )
  ),
  K = kStart * (nAgent/0.9) * 3,
  P = 1,
  FF = Cobb_Douglas$new(alpha = 0.36),
  PSI = NoShock$new(mu = 1),
  XI  = NoShock$new(mu = 1),
  delta = 0.025
)
testKs3$updatePolicies(k = testKs3$K/testKs3$L, law_k = function(k) k, tol = 0.04)
testKs3$plotDist()
testKs3$Agents[[1]]$plotPolicy()
testKs3$Agents[[2]]$plotPolicy()
testKs3$findWealthSS(verbose = T)
testKs3$optimizeStohastic(
  tol = .02, 
  tol_policy = 0.04,
  tol_ss = .02
)

#----------------------------------------------------------------------

testKs7 <- KS_Economy$new(
  Agents = list(
    # Type 1
    AgentType$new(
      M = 15,
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
      beta = 0.99 # discount factor
    ),
    # Type 2
    AgentType$new(
      M = 15,
      p = 1, # personal productivity
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
      beta = 0.986 # discount factor
    ),
    # Type 3
    AgentType$new(
      M = 15,
      p = 1, # personal productivity
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
      beta = 0.983 # discount factor
    ),
    # Type 4
    AgentType$new(
      M = 15,
      p = 1, # personal productivity
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
      beta = 0.980 # discount factor
    ),
    # Type 5
    AgentType$new(
      M = 15,
      p = 1, # personal productivity
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
      beta = 0.976 # discount factor
    ),
    # Type 6
    AgentType$new(
      M = 15,
      p = 1, # personal productivity
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
      beta = 0.973 # discount factor
    ),
    # Type 7
    AgentType$new(
      M = 15,
      p = 1, # personal productivity
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
      beta = 0.970 # discount factor
    )
  ), # list of agents
  P = 1, # aggregate productivity multiplier
  FF = Cobb_Douglas$new(alpha = 0.36), # production function
  PSI = NoShock$new(mu = 1), # aggregate permanent productivity schock
  XI  = NoShock$new(mu = 1), # aggregate transitory productivity schock
  delta = 0.025, # depreciation
  Z = 1 # TFP
)
testKs7$optimizeStohastic(
  tol = .01, 
  tol_policy = 0.01,
  tol_ss = .02
)
