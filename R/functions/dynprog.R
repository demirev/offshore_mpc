library(splines)
library(gam)

discretize_m <- function(
  P_bounds = c(0.5,1.5),
  p_bounds = c(0.5,1.5),
  alpha = 0.3, 
  beta = 0.98, 
  delta = 0.4, 
  prodf = Cobb_Douglas$new(alpha = alpha),
  above = 1.2,
  D = 0.00625,
  num_out = 500
) {
  highest_k <- ((alpha*beta)/(1 - beta*(1 - delta)))^(1/(1 - alpha))/P_bounds[1]
  max_r <- prodf$MPKk(highest_k) 
  max_m <- 1/(1 - D) * (1 - delta * max_r)
  return(seq(0, max_m, length.out = num_out))
}

discretize_k <- function(
  P_bounds = c(0.5,1.5), 
  below = 0.5, 
  above = 1.5,
  num_out = 100,
  alpha = 0.3, 
  beta = 0.98, 
  delta = 0.3
) {
  SS_k <- ((alpha*beta)/(1 - beta*(1 - delta)))^(1/(1 - alpha))
  SS_k_high <- (SS_k / P_bounds[1])*above # for lowest productivity
  SS_k_low  <- (SS_k / P_bounds[2])*below # for highest productivity
  
  return(seq(SS_k_low, SS_k_high, length.out = num_out))
}

fit_loess <- function(V) {
  if (length(unique(V$k)) > 1) {
    fit <- loess(
      action ~ m + k, 
      data = V, control = loess.control(surface = "direct")
    )
  } else {
    fit <- loess(
      action ~ m, 
      data = V, control = loess.control(surface = "direct")
    )
  }
  
  pred_func <- function(m,k){
    if (length(m) > 1 & length(k) == 1) {
      k <- rep(k, length(m))
    }
    predict(fit, data.frame(m, k, value = 0, action = 0))
  }
  return(pred_func)
}

fit_spline <- function(V, df = 3) {
  if (length(unique(V$k)) > 1) {
    fit <- gam(
      action ~ s(m, df = df) + s(k, df = df), 
      data = V
    )
  } else {
    fit <- gam(
      action ~  s(m, df = 10), 
      data = V
    )
  }
  
  pred_func <- function(m,k){
    if (length(m) > 1 & length(k) == 1) {
      k <- rep(k, length(m))
    }
    predict(fit, data.frame(m, k, value = 0, action = 0))
  }
  return(pred_func)
}

fit_loess2 <- function(V, degree = 1) {
  if (length(unique(V$k)) > 1) {
    fit <- gam(
      action ~ lo(m, degree = degree) + lo(k, degree = degree), 
      data = V
    )
  } else {
    fit <- gam(
      action ~ lo(m, degree = 2), 
      data = V
    )
  }
  
  pred_func <- function(m,k){
    if (length(m) > 1 & length(k) == 1) {
      k <- rep(k, length(m))
    }
    predict(fit, data.frame(m, k, value = 0, action = 0))
  }
  return(pred_func)
}

opt_pol <- function(
  m,
  k,
  V,
  delta = 0.3,
  D = 0.00625,
  beta = 0.985,
  polf,
  k_law = function(k) {return(k)},
  prodf = Cobb_Douglas$new(alpha = 0.3),
  utilf = Iso_Elastic$new(rho = 1),
  psi = LogNormalShock$new(),
  xi  = LogNormalShock$new(),
  cgrid = 40,
  ndraw = 10
) {
  
  #polf  <- fit_policy(V)
  #browser()
  # action-space grid
  c_seq <- seq(0, m, length.out = cgrid) # cannot consume more than m
  
  # law of motion
  future_k <- k_law(k)
  r <- prodf$MPKk(future_k)
  
  # stohastic grid
  psi_r <- psi$draw(ndraw)
  colnames(psi_r) <- c("psi","psi_p")
  xi_r  <- xi$draw(ndraw)
  colnames(xi_r)  <- c("xi", "xi_p")
  stoh_grid <- psi_r %>%
    crossing(xi_r)
  stoh_grid$prob <- stoh_grid$psi_p * stoh_grid$xi_p
  
  foc <- function(c) {
    # future assets for all psi-xi realizations
    future_m <- (1 - delta + r)*(m - c)/((1 - D)*stoh_grid$psi) + stoh_grid$xi
    future_m[future_m < 0] <- 0 # negative assets don't make sense
    # query policy at given future value
    future_c <- polf(future_m, future_k)
    # take expectation
    future_c <- sum(future_c * stoh_grid$prob)
    # for numerical stability
    future_c[future_c < 0] <- 1e-16 # small but positive
    # expected marginal utility
    exp_mu <- utilf$MU(future_c)
    # Euler equation
    return(utilf$MU(c) - beta*(1 - delta + r)*exp_mu)
  }
  
  if (foc(0) < 0) {
    return(0) # corner - consume zero
  } else if (foc(c_seq[length(c_seq)]) > 0) {
    return(c_seq[length(c_seq)]) # corner - consume all
  } else {
    opt_cons <- uniroot(foc, c(0,c_seq[length(c_seq)]))$root # Brent
    return(opt_cons)
  }
  
}

pf_iter <- function(
  # iteration parameters
  tol = 2e-3, 
  maxiter = 50,
  fit_policy = fit_spline,
  verbose = T,
  # model parameters
  alpha = 0.36,
  beta  = 0.985,
  delta = 0.025,
  D = 0.00625,
  # model functions
  k_law = function(k) {return(k)},
  prodf = Cobb_Douglas$new(alpha = 0.36),
  utilf = Iso_Elastic$new(rho = 1),
  psi = NoShock$new(const = 1),
  xi  = NoShock$new(const = 0),
  #psi = LogNormalShock$new(),
  #xi  = NormalShock$new(),
  # optimal action choice parameters
  cgrid = 40,
  ndraw = 500,
  # discretization parameters
  P_bounds_m = c(0.5, 1.5),
  p_bounds_m = c(0.5, 1.5),
  num_out_m = 80,
  above_m = 1.2,
  P_bounds_k = c(1,1),
  below_k = 1,
  above_k = 1,
  num_out_k = 1
) {
  
  # discretize state sapce
  m_seq <- discretize_m(
    P_bounds = P_bounds_m,
    p_bounds = p_bounds_m,
    alpha = alpha,
    beta = beta,
    delta = delta,
    prodf = prodf,
    above = above_m,
    D = D,
    num_out = num_out_m
  )
  k_seq <- discretize_k(
    P_bounds = P_bounds_k,
    below = below_k,
    above = above_k,
    num_out = num_out_k,
    alpha = alpha,
    beta = beta,
    delta = delta
  )
  
  # initialize optimal value table
  V <- expand.grid(m_seq, k_seq)
  names(V) <- c("m","k")
  V$value <- 0
  V$action <- V$m*0.5 # dummy initial optimal action
  
  # main iteration loop
  iter <- 0
  diff <- Inf
  while (diff > tol & iter <= maxiter) {
    
    # recalculate policy
    oldact <- V$action
    policy <- fit_policy(V)
    
    # optimal action with current continuation policy
    newacct <- sapply(
      seq(nrow(V)), # for each state
      function(i) {
        opt <- opt_pol(
          m = V$m[i], 
          k = V$k[i], 
          V = V,
          delta = delta,
          D = D,
          beta = beta,
          polf = policy,
          k_law = k_law,
          prodf = prodf,
          utilf = utilf,
          psi = psi, 
          xi = xi,
          cgrid = cgrid,
          ndraw = ndraw
        ) # calculate new optimal action, given current policy at next state
        V$action[i] <<- opt 
        return(opt)
      })
    
    # check convergence
    diff <- max(abs(V$action - oldact))
    iter <- iter + 1
    if (verbose) cat("Completed iteration", iter, "Diff:", diff, "\n")
  }
  
  return(list(VTable = V, policy = fit_policy(V)))
}


policy_plot <- function(
  policy, 
  state_space,
  main = "Policy Plot",
  ylb = "c/m",
  xlb = "m"
) {
  
  state_space$action <- policy(state_space$m, state_space$k)
  
  if (length(unique(state_space$k)) == 1) {
    p <- ggplot(data = state_space, aes(x = m, y = action/m)) +
      geom_path() +
      theme_bw() + 
      ggtitle(main) +
      xlab(xlb) + 
      ylab(ylb)
  } else {
    stop("You haven't written this yet")
  }
  p 
}

# test1 <- pf_iter() # default settings - no shocks
# 
# test2 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.01), 
#   psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 200
# )
# test2_2 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.1), 
#   psi = LogNormalShock$new(sigma = 0.1), 
#   ndraw = 200
# )
# 
# test3 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.01), 
#   #psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 500
# )
# test3_2 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.01), 
#   #psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 1500
# )
# test3_3 <- pf_iter(
#   xi = NormalShock$new(sigma = 0.1), 
#   #psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 500
# )
# 
# test4 <- pf_iter(
#   #xi = NormalShock$new(sigma = 0.01), 
#   psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 500
# )
# test4_2 <- pf_iter(
#   #xi = NormalShock$new(sigma = 0.01), 
#   psi = LogNormalShock$new(sigma = 0.01), 
#   ndraw = 1500
# )
# test4_3 <- pf_iter(
#   #xi = NormalShock$new(sigma = 0.01), 
#   psi = LogNormalShock$new(sigma = 0.1), 
#   ndraw = 500
# )

