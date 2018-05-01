library(R6)
library(tidyverse)

NoShock <- R6Class(
  "Same methods as logNormal but returns const always",
  public = list(
    const = 1,
    
    initialize = function(const = 1) {
      self$const <- const
    },
    mean = function() {
      return(self$const)
    },
    var = function() {
      0
    },
    draw = function(n) {
      return(tibble(value = self$const, prob = 1))
    }
  )
)

NormalShock <- R6Class(
  "Normal Schock Process",
  public = list(
    mu = 0,
    sigma = 1,
    
    initialize = function(mu = 0, sigma = 1) {
      self$mu <- mu
      self$sigma <- sigma
    },
    mean = function() {
      return(self$mu)
    },
    var = function() {
      return(self$sigima)
    },
    draw = function(n) {
      # draw n realizations
      value <- sort(rnorm(n, mean = self$mu, sd = self$sigma))
      # cdf for each
      probs <- pnorm(value, mean = self$mu, sd = self$sigma)
      # will be added later to sum to 1
      bonus <- 1 - probs[length(probs)]
      # discretize
      probs <- probs - c(0, probs[-length(probs)]) # banded probability
      probs[length(probs)] <- probs[length(probs)] + bonus 
      return(tibble(value = value, prob = probs))
    }
  )
)

LogNormalShock <- R6Class(
  "Log Normal Schock Process",
  public = list(
    mu = 0,
    sigma = 1,
    
    initialize = function(mu = 0, sigma = 1) {
      self$mu <- mu
      self$sigma <- sigma
    },
    mean = function() {
      exp(self$mu + self$sigma^2/2)
    },
    var = function() {
      (exp(self$sigma^2 - 1))*exp(2*self$mu + self$sigma^2)
    },
    draw = function(n) {
      # draw n realizations
      value <- sort(rlnorm(n, meanlog = self$mu, sdlog = self$sigma))
      # cdf for each
      probs <- plnorm(value, meanlog = self$mu, sdlog = self$sigma)
      # will be added later to sum to 1
      bonus <- 1 - probs[length(probs)]
      # discretize
      probs <- probs - c(0, probs[-length(probs)]) # banded probability
      probs[length(probs)] <- probs[length(probs)] + bonus 
      return(tibble(value = value, prob = probs))
    }
  )
)

Cobb_Douglas <- R6Class(
  "Cobb-Douglas Production function",
  public = list(
    alpha = 0.3,
    initialize = function(alpha) {
      self$alpha = alpha
    },
    FF = function(K, L, l = 1, Z = 1) {
      return(Z*K^(self$alpha)*l*L^(1 - self$alpha))
    },
    MPK = function(K, L, l = 1, Z = 1) {
      return(self$alpha * Z * (K/(l*L))^(self$alpha - 1))
    },
    MPL = function(K, L, l = 1, Z = 1) {
      return((1 - self$alpha) * Z * (K/(l*L))^(self$alpha))
    },
    MPKk = function(k, Z = 1) {
      return(self$alpha * Z * (k)^(self$alpha - 1))
    },
    MPLk = function(k, Z = 1) {
      return((1 - self$alpha) * Z * (k)^(self$alpha))
    }
  )
)

Iso_Elastic <- R6Class(
  "Iso-Elastic Utility",
  public = list(
    rho = 1,
    initialize = function(rho) {
      self$rho = rho
    },
    U = function(C) {
      if (self$rho == 1) {
        return(log(C))
      } else {
        return((C^(1 - self$rho) - 1)/(1 - self$rho))
      }
    },
    MU = function(C) {
      return(1/(C^self$rho))
    }
  )
)

KS_Economy <- R6Class(
  "Krussell - Smith Model of the Economy",
  public = list(
    agents = NULL, # list of agents
    K = 0, # aggregate kapital
    P = 1, # aggregate productivity multiplier
    l = sum(lapply(agents, function(x) x$l)), # labor supply
    FF = function(K,L,Z,alpha) {Z*K^alpha*L^(1 - alpha)}, # production function
    PSI = function() {return(1)}, # aggregate permanent productivity schock
    XI  = function() {return(1)}, # aggregate transitory productivity schock
    r = NULL, # return to capital
    W = NULL, # wage rate
    
    initialize = function(K, P, l, FF, PSI, XI) {
      # initial values of every field
      self$dimensions <- lout
    },
    
    step = function() {
      # simulate a single period
      
      # 1. Produce
      
      # 2. Invest
    },
    
    calc_r = function() {
      
    },
    
    calc_W = function() {
      
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

