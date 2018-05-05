library(R6)
library(tidyverse)

KS_Economy <- R6Class(
  "Krussell - Smith Model of the Economy",
  public = list(
    Agents = NULL, # list of agents
    K = 0, # aggregate kapital
    P = 1, # aggregate productivity multiplier
    FF = Cobb_Douglas$new(alpha = 0.36), # production function
    PSI = NoShock$new(mu = 1), # aggregate permanent productivity schock
    XI  = NoShock$new(mu = 1), # aggregate transitory productivity schock
    r = NULL, # return to capital
    W = NULL, # wage rate
    L = NULL, # total labor units
    Z = NULL, # TFP
    delta = NULL, # depreciation
    nagents = NULL, # number of agents being simulated
    
    initialize = function(
      Agents = NULL, # list of agents
      K = 0, # aggregate kapital
      P = 1, # aggregate productivity multiplier
      #l = 1, # labor supply
      FF = Cobb_Douglas$new(alpha = 0.36), # production function
      PSI = NoShock$new(mu = 1), # aggregate permanent productivity schock
      XI  = NoShock$new(mu = 1), # aggregate transitory productivity schock
      delta = 0.025, # depreciation
      Z = 1 # TFP
    ) {
      # initial values of every field
      self$Agents <- Agents
      self$K <- K
      self$P <- P
      self$FF <- FF
      self$PSI <- PSI
      self$Z <- Z
      self$delta <- delta
      
      self$L <- self$update_L()
      self$nagents <- sum(sapply(Agents, function(agent) agent$n))
    },
    
    update_L = function() {
      # calculates total labor amount
      L <- self$Agents %>%
        sapply(
          function(agent) {
            sum(agent$wallets$L * agent$wallets$l)
          }
        ) %>%
        sum
      return(L)
    },
    
    quantileSummary = function(probs = seq(0,1,0.1), descending = T) {
      allWealth <- self$Agents %>%
        lapply(
          function(agent) {
            agent$wallets$M
          }
        ) %>%
        reduce(c) %>%
        sort(decreasing = descending) 
      
      cumWealth <- allWealth %>%
        cumsum %>%
        quantile(probs = probs)
      
      return(cumWealth / sum(allWealth))
    },
    
    step = function(verbose = T, predraw_xi = 30, predraw_psi = 30) {
      # simulate a single period
      #browser()
      # 1. Produce
      self$P  <- self$P * self$PSI$draw(1)$value # aggregate permanent shock
      trShock <- self$XI$draw(1)$value # aggreagete transitory shock
      # factor prices
      self$W = self$FF$MPL(
        K = self$K, 
        L = self$L * self$P * trShock, 
        l = 1, 
        Z = self$Z
      )
      self$r = self$FF$MPK(
        K = self$K, 
        L = self$L * self$P * trShock, 
        l = 1, 
        Z = self$Z
      )

      # agents receive income
      distribute_income <- self$Agents %>%
        lapply(
          function(agent) {
            agent$receive(
              W = self$W, 
              R = (1 - self$delta + self$r)
            )
          }
        )
      
      # 2. Invest and choose labor
      
      # next period aggregate capital is the sum of all investments
      self$K <- self$Agents %>%
        lapply(
          function(agent) {
            agent$invest(k = self$K / self$nagents)
          }
        ) %>%
        reduce(c) %>%
        sum
      
      # agents can potentially also change their income
      laborDecision <- self$Agents %>%
        lapply(
          function(agent) {
            agent$choose_l() # a dummy function - not part of the model for now
          }
        )
      self$L <- self$update_L()
      
      # 3. Death
      # sum the assets of the dead...
      bequests <- self$Agents %>%
        lapply(
          function(agent) {
            agent$die()
          }
        ) %>%
        reduce(c) %>%
        sum
      # ...and distribute them across survivors (and newborns)
      heirs <- self$Agents %>%
        sapply(
          function(agent) {
            agent$inherit(bequests / self$nagents) # evenly distributed
          }
        )
      
      if (verbose) {
        #cat("Period Simulated \n")
        print(self$quantileSummary())
      }
      return(T)
    },
    
    wealth_ss = function(tol = 0.0002, maxit = 6000, verbose = T) {
      wealthdist <- self$quantileSummary() # initialize wealth distr
      diff <- 2*tol # initialize higher than tol
      it <- 0
      while (diff >= tol & it <= maxit) {
        # simulate 1 period
        S <- self$step(verbose = F)
        # calculate wealth distribution
        newdist <- self$quantileSummary()
        # track change
        diff <- sum(abs(wealthdist - newdist))
        wealthdist <- newdist
        # print progress
        if (verbose & it %% 100 == 0) {
          cat("Current Difference:", round(diff,4), "\n")
          print(wealthdist)
        }
        it <- it + 1
      }
      return(TRUE)
    }
  )
)

AgentType <- R6Class(
  "Agent in a Krussell - Smith Economy",
  public = list(
    psi = LogNormalShock$new(sigma = 0.025), # productivity permanent shock
    xi = EmploymentShock$new(sigma = 0.04), # productivity transitory shock
    policy = NULL, #function(m,k) return(0.5*m), # consumption policy
    D = 0.00625, # probability of death
    beta = 0.985, # discount factor
    n = 1000, # number of agents for simulation
    wallets = NULL, # holds individual information for simulating
    
    initialize = function(
      M = NULL, # resources on hand not-nromalized
      p = 1, # personal productivity
      l = 1, # fraction of labor supplied
      L = 1, # units of labor available
      psi = LogNormalShock$new(sigma = 0.025), # productivity permanent shock
      xi = EmploymentShock$new(sigma = 0.04), # productivity transitory shock
      pol = function(m,k) return(0.5*m),
      D = 0.00625, # probability of death
      beta = 0.985, # discount factor
      n = 1000
    ) {
      # assign fields
      self$psi <- psi
      self$xi  <- xi
      self$policy <- pol
      self$D <- D
      self$beta <- beta
      self$n
      
      self$wallets <- data.frame(
        M = rep(M, n),
        p = rep(p, n),
        pW = rep(p, n),
        l = rep(l, n),
        L = rep(L, n)
      )
    },
    
    receive = function(W, R) {
      # consumption policy
      xi_val <- self$xi$sample(self$n) # to do - tau, l
      self$wallets$M <- self$wallets$M * R + self$wallets$p * xi_val * W
      self$wallets$pW <- self$wallets$p * W
      return(self$wallets$M)
    },
    
    invest = function(k) {
      psi_val <- self$psi$sample(self$n)
      m <- self$wallets$M / (self$wallets$pW) # normalize resources
      c <- self$policy(m = m, k = rep(k,self$n)) # choose consumption
      self$wallets$M <- self$wallets$M - c * self$wallets$pW # reduce resources
      self$wallets$p <- self$wallets$p * psi_val # persistent productivity shock
      return(self$wallets$M)
    },
    
    choose_l = function(k) {
      return(self$wallets$l * self$wallets$L) # dummy method
      # maybe intriduced in the future
    },
    
    die = function() {
      dead <- (runif(self$n) < self$D)
      bequest <- self$wallets$M * dead
      self$wallets$M[dead] <- 0
      self$wallets$p[dead] <- 1
      self$wallets$pW[dead] <- 1
      
      return(bequest)
    },
    
    inherit = function(bequest) {
      # increase resources by inheritence
      self$wallets$M <- self$wallets$M + bequest
      return(T)
    }
  )
)

testKs <- KS_Economy$new(
  Agents = list(
    AgentType$new(
      M = 10.7, 
      p = 1, 
      l = 1/0.9,
      L = 1,
      psi = LogNormalShock$new(sigma = 0.025), 
      xi = EmploymentShock$new(sigma = 0.04), 
      pol = test6_3$policy,
      D = 0.00625,
      beta = 0.985,
      n = 10000
    )
  ),
  K = 10.7 * 1000, 
  P = 1, 
  FF = Cobb_Douglas$new(alpha = 0.36),
  PSI = NoShock$new(mu = 1),
  XI  = NoShock$new(mu = 1),
  delta = 0.025 
)
