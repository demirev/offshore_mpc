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
    
    quantileSummary = function(probs = seq(0,1,0.1), 
                               descending = T, verbose = F) {
      allWealth <- self$Agents %>%
        lapply(
          function(agent) {
            agent$wallets$M / agent$wallets$pW # permanent income
          }
        ) %>%
        reduce(c) %>%
        sort(decreasing = descending) 
      
      cumWealth <- allWealth %>%
        cumsum %>%
        quantile(probs = probs)
      
      if (verbose) hist(allWealth)
      
      return(
        list(
          quantiles = cumWealth / sum(allWealth),
          x = allWealth
        )
      )
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
        print(self$quantileSummary()$quantiles)
      }
      
      return(self$K / self$L)
    },
    
    findWealthSS = function(
      tol = 0.02, 
      maxit = 6000, 
      minit = 500, 
      burnin = 200, 
      verbose = T
    ) {
      # simulates the economy until a steady wealth distribution is reached
      
      wealthdist <- self$quantileSummary()$quantiles # initialize wealth distr
      wealthx    <- self$quantileSummary()$x # test
      
      k_history  <- vector("numeric", maxit + 1)
      diff <- 2*tol # initialize higher than tol
      it <- 0
      while ((diff >= tol & it <= maxit) | (it <= minit)) {
        
        # simulate 1 period
        k_history[it + 1] <- self$step(verbose = F, fixK = F)
        
        # print progress
        if (it %% 100 == 0) {
          # calculate wealth distribution
          newdist <- self$quantileSummary()$quantiles
          newx    <- self$quantileSummary()$x
          
          # track change
          diff <- sum(abs(wealthdist - newdist))
          wealthdist <- newdist
          kstest <- ks.test(wealthx, newx)
          wealthx <- newx
          
          # progress message
          if (verbose) {
            cat("Current Difference:", round(diff,4), "\n")
            print(self$quantileSummary(verbose = T)$quantiles)
            print(kstest)
            KoL <- self$K / self$L
            cat("Current K/L:", KoL, "\n")
          }
          
        }
        
        it <- it + 1
      }
      return(k_history[(burnin - 1):it]) # returns everything but the burnin
    },
    
    updatePolicies = function(
      k, 
      law_k, 
      verbose = T, 
      fit_policy = fit_spline,
      onlySS = F,
      max_m = NULL,
      num_out = NULL,
      tol = 4e-2,
      firstpass = TRUE
    ) {
      #browser()
      for (agent_no in seq(length(self$Agents))) {
        if (firstpass) {
          start_from <- max(1,agent_no - 1)
        } else {
          start_from <- agent_no
        }
        
        self$Agents[[agent_no]]$updatePolicy(
          alpha = self$FF$alpha,
          delta = self$delta, 
          prodf = self$FF, 
          ss_k = k,
          start_action = self$Agents[[start_from]]$QTable$action, # hot start value iteration
          law_k = law_k, 
          fit_policy = fit_policy, 
          onlySS = onlySS, 
          max_m = ifelse(is.null(max_m), self$Agents[[agent_no]]$max_m, max_m),
          num_out = ifelse(is.null(num_out), self$Agents[[agent_no]]$num_out, num_out),
          tol = tol
        )
      }
      # new_policies <- self$Agents %>%
      #   lapply(
      #     function(agent) {
      #       agent$updatePolicy(
      #         alpha = self$FF$alpha,
      #         delta = self$delta, 
      #         prodf = self$FF, 
      #         ss_k = k,
      #         start_action = self$Agents[[1]]$QTable$action, # hot start value iteration
      #         law_k = law_k, 
      #         fit_policy = fit_policy, 
      #         onlySS = onlySS, 
      #         max_m = ifelse(is.null(max_m), agent$max_m, max_m),
      #         num_out = ifelse(is.null(num_out), agent$num_out, num_out),
      #         tol = tol
      #       )
      #     }
      #   )
      
      if (verbose) {
        for (agent in self$Agents) {
          p <- agent$plotPolicy()
          print(p)
          #print(agent$QTable)
        }
      }
      
      return(T)
    },
    
    optimizeStohastic = function(
      k = 38.9, 
      tol = .02, 
      tol_policy = 4e-2,
      tol_ss = .02,
      lrate = .3, 
      max_m = NULL,
      num_out = NULL,
      fit_policy = fit_spline,
      verbose = T,
      probs = seq(0,1,0.1)
    ) {
      # iterates between policy estimation and simulation until estimated
      # law of motion matches the actual one
      converged <- F
      
      # initialzie law of motion guess
      a0 <- 0
      a1 <- 0.95
      law_k <- function(k, b0 = a0, b1 = a1) {
        return(exp(b0 + b1*log(k)))
      }
      self$K <- k * self$L
      
      #browser()
      firstpass <- TRUE # at first pass the first agent's policy will be cold
      # started, and every other agent will start from that policy
      
      while (!converged) {
        
        # 0. Reset m
        newWealth <- self$Agents %>%
          lapply(
            function(agent) {
              agent$wallets$M <- 0#k
            }
          )
        
        # 1. Calculate Policy Given steady state guess
        update_policies <- self$updatePolicies(
          k = k, 
          law_k = law_k, 
          verbose = T, 
          fit_policy = fit_policy, 
          onlySS = F, 
          max_m = max_m, 
          num_out = num_out, 
          tol = tol_policy,
          firstpass = firstpass
        )
        firstpass <- FALSE # after the first pass each agent's policy will be
        # hot started from their pervious policy
        
        # 2. Simulate
        wealth_ss <- self$findWealthSS(verbose = T, minit = 1000, tol = tol_ss)
        
        # 3. update law of motion guess
        dep <- log(wealth_ss[2:length(wealth_ss)]) # everything but the first one
        lag <- log(wealth_ss[1:(length(wealth_ss) - 1)]) # fit against 1-lag
        reg <- lm(dep ~ lag)
        reg$coefficients[is.na(reg$coefficients)] <- 0 
        
        diff <- max(
          abs(a0 - reg$coefficients[1]), 
          abs(a1 - reg$ceofficients[2])
        ) # update magnitude
        
        if (diff < tol) {
          converged <- TRUE
        } else {
          # update
          a0 <- (1 - lrate)*a0 + lrate * reg$coefficients[1]
          a1 <- (1 - lrate)*a1 + lrate * reg$coefficients[2]
          law_k <- function(k, b0 = a0, b1 = a1) {
            return(exp(b0 + b1*log(k)))
          }
        }
        
        self$K <- mean(wealth_ss) * self$L # for next iteration
        
        # 4. Message
        if (verbose) {
          cat(
            "Iteration complete. Diff:", diff, 
            "coefficients:", a0, a1, 
            "\n"
          )
        }
      }
      
      return(self$quantileSummary(probs = probs)$quantiles)
      
    },
    
    calcMPC = function(
      fit_policy = fit_spline,
      at_k = self$K / self$L,
      eps = 0.02 # epsilon increase in m leads to x% increase in spending
    ) {
      mpc_list <- self$Agents %>%
        lapply(
          function(agent) {
            QTable <- agent$QTable %>%
              filter(k == at_k) # only the selected capital value
            
            if (nrow(QTable) == 0) {
              stop("This value of capital is not")
            }
            
            policy <- fit_policy(
              agent$QTable, 
              df_m = length(unique(agent$QTable$m))
            ) # interpolate again using many degrees of freedom
            
            m <- agent$wallets$M / agent$wallets$pW
            
            mpc_vector <- (
              policy(m + eps, rep(at_k, nrow(agent$wallets))) - 
              policy(m, rep(at_k, nrow(agent$wallets)))
            )/eps
            
            return(mpc_vector)
          }
        )
      
      mpc <- mpc_list %>%
        reduce(c) %>%
        mean
      
      return(list(mpc = mpc, mpc_list = mpc_list))
    },
    
    plotDist = function(nbin = 30, padTo = 120, 
                        type = "plotly", ylimC = c(0,8), ylimW = c(0,0.25)) {
      
      wealth <- self$Agents %>%
        lapply(
          function(agent) {
            tibble(
              m = agent$wallets$M / agent$wallets$pW,
              beta = agent$beta
            )
          }
        ) %>%
        reduce(rbind)
      
      policies <- self$Agents %>%
        lapply(
          function(agent) {
            maxm <- max(agent$QTable$m)
            mpad <- seq(maxm, padTo, 1)[-1]
            kpad <- rep(agent$QTable$k[1], length(mpad))
            tibble(
              m = c(agent$QTable$m, mpad),
              c = c(agent$QTable$action, agent$policy(m = mpad, k = kpad)),
              beta = agent$beta
            )
          }
        ) %>%
        reduce(rbind)
      
      if (type == "plotly") {
        p <- plot_ly() %>%
          add_trace(
            x = policies$m, y = policies$c, type = 'scatter', 
            mode = 'lines', name = 'consumption',  yaxis = 'y2',
            line = list(color = '#45171D'), 
            hoverinfo = "text",
            text = policies$c
          ) %>%
          add_trace(
            x = wealth$m, type = 'histogram', name = 'wealth', nbinsx = nbin,
            marker = list(color = '#C9EFF9'),
            histnorm = "probability"
          ) %>%
          layout(
            title = 'Consumption and Wealth',
            xaxis = list(
              title = "m"#,
              #range = c(0, 60)
            ),
            yaxis2 = list(
              side = 'right', 
              title = 'c', 
              overlaying = "y",
              showgrid = FALSE, 
              zeroline = FALSE,
              range = ylimC
            ),
            yaxis = list(
              side = 'left', 
              title = '', 
              showgrid = FALSE, 
              zeroline = FALSE,
              range = ylimW
            )
          )
      } else if (type == "ggplot") {
        policies$panel = "Policy"
        wealth$panel = "Wealth"
        p <- ggplot() + 
          facet_grid(panel ~ ., scales = "free") + 
          geom_line(data = policies, aes(x = m, y = c), col = "red") +
          geom_histogram(data = wealth, aes(x = m), alpha = .5, bins = nbin) +
          theme_bw()
      }
      return(p)
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
    wallets = NULL, # holds individual information for simulating,
    max_m = NULL, # maximal m value for polict iteration (extrapolation after that)
    num_out = NULL, # number of grid points on the m-grid
    
    initialize = function(
      M = 15, # resources on hand not-nromalized
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
      max_m = 35, # maximal m value for polict iteration (extrapolation after that)
      num_out = 45, # number of grid points on the m-grid
      n = 1000
    ) {
      
      if (M <= 0) stop("Please initialize agents with positive wealth")
      
      # assign fields
      self$psi <- psi
      self$xi  <- xi
      self$utilf <- utilf
      self$policy <- pol
      self$QTable <- QTable
      self$D <- D
      self$beta <- beta
      self$n <- n
      self$max_m <- max_m
      self$num_out <- num_out
      
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
    
    updatePolicy = function(
      alpha, 
      delta, 
      prodf, 
      ss_k, 
      start_action = NULL,
      law_k,
      fit_policy = fit_spline,
      onlySS = F,
      max_m = self$max_m,
      num_out = self$num_out,
      tol = 4e-2
    ) {
      # given parameters of the economy solves for optimal policy
      # this version assumes no aggregate shocks, so only parameter is 
      # steady state capital
      #browser()
      if (!onlySS) {
        k_seq <- c(0.2*ss_k, 0.6*ss_k, ss_k, 1.5*ss_k, 3*ss_k)
      } else {
        k_seq <- ss_k
      }
      
      newPol <- pf_iter(
        # iteration parameters
        tol = 4e-2, 
        maxiter = 50,
        fit_policy = fit_policy, 
        verbose = T, 
        action = start_action,
        # model parameters
        alpha = alpha,
        beta  = self$beta,
        delta = delta,
        D = self$D,
        # model functions
        k_law = law_k, #function(k) {return(k)}, # always in steady state
        prodf = prodf,
        utilf = self$utilf,
        psi = self$psi,
        xi  = self$xi,
        # optimal action choice parameters
        cgrid = 40,
        ndraw = 20,
        # discretization parameters
        m_seq = discretize_m(
          max_m = max_m,
          num_out = num_out
        ),
        k_seq = k_seq # around s.s. k
      )
      
      self$policy <- newPol$policy
      self$QTable <- newPol$QTable
    },
    
    plotPolicy = function() {
      p <- ggplot() + 
        geom_line(
          data = self$QTable, 
          aes(x = m, y = action, group = k, color = k)
        ) +
        labs(title = expression(paste(beta, "=", self$beta))) +
        theme_bw() +
        scale_colour_gradient(low = "gray", high = "black")
      return(p)
    }
  )
)
