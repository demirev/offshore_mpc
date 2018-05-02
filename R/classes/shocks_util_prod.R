library(R6)
library(tidyverse)

NoShock <- R6Class(
  "Same methods as logNormal but returns const always",
  public = list(
    mu = 1,
    sigma = 0, # its a constant
    
    initialize = function(mu = 1, sigma = 0) {
      self$mu <- mu
    },
    mean = function() {
      return(self$const)
    },
    var = function() {
      return(self$sigma^2)
    },
    draw = function(n) {
      return(tibble(value = self$mu, prob = 1))
    }
  )
)

NormalShock <- R6Class(
  "Normal Schock Process",
  inherit = NoShock,
  public = list(
    mu = 0,
    sigma = 1,
    initialize = function(mu = 0, sigma = 1) {
      self$mu <- mu
      self$sigma <- sigma
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
  inherit = NormalShock,
  public = list(
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

EmploymentShock <- R6Class(
  "Employment per the Caroll specification",
  public = list(
    mu_ins = 0.15, # insurance
    #tau = 0, # tax
    #l = 1, # portion of labor supplied
    OM = 0.07, # probability of unemployment
    theta = NULL, #
    
    initialize = function(
      mu = 0, 
      sigma = 1,
      mu_ins = 0.15,
      #tau = 3,
      #l = 1,
      OM = 0.07
    ) {
      self$theta <- LogNormalShock$new(mu, sigma)
      self$mu_ins <- mu_ins
      #self$tau = tau
      #self$l = l
      self$OM = OM
    },
    mean = function(l = 1, tau = 0.02) {
      return(self$mu_ins*self$OM + (1 - self$OM)*self$theta$mean()*l*(1 - tau))
    },
    var = function() {
      return(self$theta$var())
    },
    draw = function(n, l = 1, tau = 0.02) {
      # For employment:
      theta_draw <- self$theta$draw(n)
      value <- theta_draw$value * l * (1 - tau)
      probs <- theta_draw$prob * (1 - self$OM)
      # For unemployment
      value <- c(self$mu_ins, value)
      probs <- c(self$OM, probs)
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