library(R6)
library(tidyverse)

KS_Economy <- R6Class(
  "Krussell - Smith Model of the Economy",
  public = list(
    agents = NULL, # list of agents
    K = 0, # aggregate kapital
    P = 1, # aggregate productivity multiplier
    FF = Cobb_Douglas$new(), # production function
    PSI = NoShock$new(mu = 1), # aggregate permanent productivity schock
    XI  = NoShock$new(mu = 1), # aggregate transitory productivity schock
    r = NULL, # return to capital
    W = NULL, # wage rate
    L = NULL, # total labor units
    
    initialize = function(
      agents = NULL, # list of agents
      K = 0, # aggregate kapital
      P = 1, # aggregate productivity multiplier
      l = 1, # labor supply
      FF = Cobb_Douglas$new(alpha = 0.36), # production function
      PSI = NoShock$new(mu = 1), # aggregate permanent productivity schock
      XI  = NoShock$new(mu = 1), # aggregate transitory productivity schock
    ) {
      # initial values of every field
      self$agents <- agents
      self$K <- K
      self$P <- P
      self$FF <- FF
      self$PSI <- PSI
      
      self$L <- sum(sapply(agents, function(x) x$L*x$l))
    },
    
    step = function() {
      # simulate a single period
      
      # 1. Produce
      self$P  <- self$P * self$PSI$draw(1)$value
      trShock <- self$XI$draw(1)$value
      W = self$FF$MPL(
        K = self$K, 
        L = self$L * self$P * trShock, 
        l = 1, 
        Z = self$Z
      )
      R = self$FF$MPK(
        K = self$K, 
        L = self$L * self$P * trShock, 
        l = 1, 
        Z = self$Z
      )
      
      distribute_income <- sapply(
        self$agents, 
        function(agent) {
          agent$receive(W, R)
        }
      )
      
      # 2. Invest
      self$K <- self$agents %>%
        sapply(
          function(agent) {
            agent$invest()
          }
        ) %>%
        sum
      
      # 3. Death
      
    }
    
    
  )
)

Agent <- R6Class(
  "Agent in a Krussell - Smith Economy",
  public = list(
    loc          = NULL, # location
    
    initialize = function() {
      self$dead        <- FALSE # it's alive
    },
    
    consume = function() {
      # consumption policy
      self$loc <- loc
      return(loc)
    }
  )
)

