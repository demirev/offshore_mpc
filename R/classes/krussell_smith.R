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
    
    step = function(verbose = T, fixK = FALSE) {
      # simulate a single period
      #browser()
      # 1. Produce
      self$P  <- self$P * self$PSI$sample(1) # aggregate permanent shock
      trShock <- self$XI$sample(1) # aggreagete transitory shock
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
      investments <- self$Agents %>%
        lapply(
          function(agent) {
            agent$invest(k = self$K / self$nagents)
          }
        ) %>%
        reduce(c) %>%
        sum
      if (!fixK) self$K <- investments
       
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
    
    wealth_ss = function(tol = 0.02, maxit = 6000, minit = 500, verbose = T) {
      # simulates the economy until a steady wealth distribution is reached
      
      wealthdist <- self$quantileSummary() # initialize wealth distr
      diff <- 2*tol # initialize higher than tol
      it <- 0
      while ((diff >= tol & it <= maxit) | (it <= minit)) {
        # simulate 1 period
        S <- self$step(verbose = F, fixK = T)
        # calculate wealth distribution
        newdist <- self$quantileSummary()
        # track change
        diff <- sum(abs(wealthdist - newdist))
        wealthdist <- newdist
        # print progress
        if (verbose & it %% 100 == 0) {
          cat("Current Difference:", round(diff,4), "\n")
          print(wealthdist)
          KoL <- self$K / self$L
          cat("Current K/L:", KoL, "\n")
        }
        it <- it + 1
      }
      return(self$K / self$L)
    },
    
    stohastic_optimization = function(
      k = 38.9, 
      tol = 0.02, 
      lrate = .3, 
      verbose = T
    ) {
      # iterates between policy estimation and simulation until the K/L implied
      # by the policy matches the one used to estimate it
      converged <- F
      
      while (!converged) {
        # 1. Calculate Policy Given steady state guess
        update_policies <- self$Agents %>%
          lapply(
            function(agent) {
              agent$updatePolicy(
                alpha = self$FF$alpha,
                delta = self$delta, 
                prodf = self$FF, 
                ss_k = k,
                start_action = agent$QTable$action # hot start value iteration
              )
            }
          )
        # 2. Simulate
        wealth_ss <- self$wealth_ss(verbose = T)
        new_k <- self$Agents %>%
          lapply(
            
            function(agent) {
              agent$wallets$M
            }
          ) %>%
          reduce(c) %>%
          sum # implied K by agents policy
        new_k <- new_k / self$L
        diff <- abs(new_k - k)
        if (diff < tol) {
          converged <- TRUE
        } else {
          if (new_k > 3*k) {
            k <- 3*k # avoid enourmous updates
          } else {
            k <- (1 - lrate) * k + lrate * new_k
          }
        }
        
        self$K <- k * self$L
        
        # 3. Message
        if (verbose) {
          cat("Iteration complete. Diff:", diff, "k's:", new_k, k, "\n")
        }
      }
      
    }
  )
)

AgentType <- R6Class(
  "Agent in a Krussell - Smith Economy",
  public = list(
    psi = LogNormalShock$new(sigma = 0.025), # productivity permanent shock
    xi = EmploymentShock$new(sigma = 0.04), # productivity transitory shock
    utilf = Iso_Elastic$new(rho = 1), # utility function
    policy = NULL, #function(m,k) return(0.5*m), # consumption policy
    QTable = data.frame(m = 0, k = 0, action = 0), # optimal action table
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
      utilf = Iso_Elastic$new(rho = 1), # utility function
      pol = function(m,k) return(0.5*m),
      QTable = data.frame(m = 0, k = 0, action = 0),
      D = 0.00625, # probability of death
      beta = 0.985, # discount factor
      n = 1000
    ) {
      # assign fields
      self$psi <- psi
      self$xi  <- xi
      self$utilf <- utilf
      self$policy <- pol
      self$QTable <- QTable
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
      # kill off (and respawn) some of the agents
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
    },
    
    updatePolicy = function(alpha, delta, prodf, ss_k, start_action = NULL) {
      # given parameters of the economy solves for optimal policy
      # this version assumes no aggregate shocks, so only parameter is 
      # steady state capital
      
      newPol <- pf_iter(
        # iteration parameters
        tol = 3e-2, 
        maxiter = 50,
        fit_policy = fit_spline,
        verbose = T,
        # model parameters
        alpha = alpha,
        beta  = self$beta,
        delta = delta,
        D = self$D,
        # model functions
        k_law = function(k) {return(k)}, # always in steady state
        prodf = prodf,
        utilf = self$utilf,
        psi = self$psi,
        xi  = self$xi,
        # optimal action choice parameters
        cgrid = 40,
        ndraw = 20,
        # discretization parameters
        m_seq = discretize_m(
          max_m = 35,
          num_out = 130
        ),
        k_seq = ss_k # steady state capital
      )
      
      self$policy <- newPol$policy
      self$QTable <- newPol$QTable
    }
  )
)

testKs <- KS_Economy$new(
  Agents = list(
    AgentType$new(
      M = 38.9, 
      p = 1, 
      l = 1/0.9,
      L = 1,
      psi = LogNormalShock$new(sigma = 0.025), 
      xi = EmploymentShock$new(sigma = 0.04), 
      utilf = Iso_Elastic$new(rho = 1),
      pol = test6_3$policy,
      QTable = test6_3$VTable,
      D = 0.00625,
      beta = 0.985,
      n = 10000
    )
  ),
  K = 38.9 * 10000, 
  P = 1, 
  FF = Cobb_Douglas$new(alpha = 0.36),
  PSI = NoShock$new(mu = 1),
  XI  = NoShock$new(mu = 1),
  delta = 0.025 
)

testKs$stohastic_optimization()
