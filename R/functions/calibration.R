lossKS <- function(target) {
  # evaluates total distance between simulated and target quantiles (generator)
  ff <- function(quantiles) {return(sum(abs(quantiles - target)))}
  return(ff)
}

lossKSQuintiles <- function(target) {
  # evaluates total distance between simulated and target quintiles (generator), 
  # assuming targets come in quintiles and the data come in deciles
  if (length(target) != 6){
    stop("lossKSQunitiles can only be used with targets that have 6 columns")
  }
  ff <- function(quantiles) {
    if (length(quantiles) != 11){
      stop("lossKSQunitiles can only be used with quantile data that have 11 columns")
    }
    # Use columns 1, 3, 5, 7, 9, 11 (0%, 20%, 40%, 60%, 80%, 100%)
    quintiles = quantiles[c(1, 3, 5, 7, 9, 11)]
    return(sum(abs(quantiles - target)))
  }
  return(ff)
}

generateKSParams <- function(
  beta_mid_span = c(0.9, 0.99), 
  beta_rng_span = c(0.01, 0.1)
) {
  # generates a single couple of beta_mid and beta_range (generator)
  param_drawer <- function() {
    beta_mid = 2
    beta_rng = 0
    # generate beta_mid and beta_rng until all betas are between 0 and 1
    while(beta_mid + beta_rng > 1 | beta_mid - beta_rng < 0) {
      beta_mid = beta_mid_span[1] + runif(1)*(beta_mid_span[2] - beta_mid_span[1])
      beta_rng = beta_rng_span[1] + runif(1)*(beta_rng_span[2] - beta_rng_span[1])
    }
    return(
      c(beta_mid = beta_mid, beta_range = beta_rng)
    )
  }
  return(param_drawer)
}

calibrate_genetic <- function(
  FUN, # function to be optimized. Returns a single object...
  lossF, # that is past to the loss function - returns a num (lower is bettter)
  individual_generator, # a function that generates a single vector of numeric hyper-parameters
  npop = 20, # size of population
  nsurvive = 6, # number of survivors per generation
  generations = 10, # number of generations to train
  tol = 1e-2, # will stop early if loss is less than this
  nparents = 3, # number of parents per children
  nchild = 10, # number of children to spawn per generation
  initial_pop = NULL, # starting population can be supplied directly
  checkpoint = NULL, # file to write results to
  recordOutput = F, # if performance is saved to file should only fitness be
  logMessages = F, # if output should be recorded
  log_file = NULL # file to write stdout to
) {
  # Performs a simple genetic optimization for hyper-parameter tuning
  
  # set up logging
  if (logMessages) {
    sink(log_file, append = T, split = T)
    # a function that closes all open sinks
    sink.reset <- function(){
      for(i in seq_len(sink.number())){
        sink(NULL)
      }
    }
    on.exit(sink.reset()) # stop sinking once we're done
  }
  
  if (nchild + nsurvive > npop - 1) stop("Reduce children or survivors")
  if (nchild > choose(nsurvive, nparents)) stop("Too many children requested")
  
  # helpers --------------------------------------------------------------------
  spawn_n <- function(n) {
    # helper - spawns n draws from individual_generaotr
    n %>%
      seq %>%
      lapply(
        function(x) {
          individual_generator()
        }
      )
  }
  
  give_birth <- function(parents) {
    # helper - combines parents to create "children"
    parents %>%
      reduce(rbind) %>%
      colMeans
  }
  
  loss_message <- function(loss, individual) {
    cat("\t Evaluated individual", individual, "Loss is", loss,"\n")
    return(loss)
  }
  
  lossWrapper <- function(obj) {
    # takes the output of FUN passes it to lossF and returns both
    loss <- lossF(obj)
    return(list(loss = loss, obj = obj))
  }
  
  checkWriter <- function(lossObj, individual) {
    # appends performance to csv (if needed)
    if (!is.null(checkpoint)) {
      writeColnames <- !file.exists(checkpoint)
      names(individual) <- paste0("param_", names(individual))
      toWrite <- c(individual, fitness = lossObj$loss)
      if (recordOutput) toWrite <- c(toWrite, lossObj$obj)
      toWrite <- toWrite %>%
        t %>% 
        as_tibble %>%
        write_csv(
          path = checkpoint, 
          append = T, 
          col_names = writeColnames
        )
    }
    return(lossObj$loss) # no need to carry lossObj$obj further
  }
  
  # generate initial population ------------------------------------------------
  if (!is.null(initial_pop)) {
    # Initial population is supplied - use it (but don't skip first evaluation)
    population <- initial_pop
    skipFirst <- F
  } else if (is.null(checkpoint)) {
    # No initial population, no checkpoint file - generate population
    population <- spawn_n(npop)
    skipFirst <- F
  } else if (!file.exists(checkpoint)) {
    # generate population but create checkpoint file
    population <- spawn_n(npop)
    #file.create(checkpoint)
    skipFirst <- F
  } else {
    # choose population from checkpoint. Don't evaluate initial
    checkpoint_file <- checkpoint %>%
      read_csv
    
    if (nrow(checkpoint_file) >= npop) {
      population <- checkpoint %>%
        read_csv %>%
        arrange(fitness) %>%
        slice(1:npop) 
      
      fitness <- population$fitness
      
      population <- population %>%
        transpose %>%
        lapply(function(member) {
          nms <- names(member)
          result <- reduce(member, c)
          names(result) <- nms
          result <- result[str_detect(nms, "param_")]
          names(result) <- str_replace(names(result), "param_", "")
          return(result)
        }) # some annoying data wrangling to get it in list of named vecotrs format
      skipFirst <- T
    } else {
      population <- spawn_n(npop)
      skipFirst <- F
    }
    
  }
  
  # main loop ------------------------------------------------------------------
  for (generation in seq(generations)) {
    print(paste0(Sys.time(), "-------------------------------------------"))
    cat("------- Training Generation", generation, "---------\n")
    print(reduce(population, rbind))
    
    if (skipFirst & generation == 1) {
      parents <- population[1:nsurvive] # if performance is loaded from 
      #checkpoint no need to reavaluate fitness on first pass
      parents_fitness <- fitness[1:nsurvive]
    } else {
      # calculate fitness
      fitness <- switch(as.character(as.numeric(generation == 1)),
        "1" = population,
        "0" = c(children, mutants) # no need to evaluate the parents twice
      ) %>%
        sapply(
          function(individual) {
            individual %>%
              FUN %>%
              lossWrapper %>%
              checkWriter(individual) %>%
              loss_message(individual)
          }
        )
      
      if (generation > 1) fitness <- c(parents_fitness, fitness)
      
      if (min(fitness) < tol) {
        cat("Converged to target. Best parameters are: \n")
        cat(population[[which.min(fitness)]])
        break()
      }
      
      # cutoff for survival
      maxfitness <- sort(fitness, decreasing = F)[nsurvive]
      
      # choose survivors      
      parents <- population[fitness <= maxfitness]
      parents_fitness <- fitness[fitness <= maxfitness]
      parents <- parents[1:nsurvive] # in case of tie
      parents_fitness <- parents_fitness[1:nsurvive]
      
      currentBest <- population[[which.min(fitness)]]
      
      cat(
        "^ Current best parameters are:", 
        currentBest, 
        "with loss:", 
        min(fitness),
        "\n"
      )
    }
        
    # determine parent 'couples'
    combos <- combn(nsurvive, nparents)
    combos <- combos[ ,sample(1:ncol(combos), nchild)]
    children <- parents %>% # start with parents
      length %>% # count them
      combn(nparents) %>% # get number of 'couples' (or triplets or whatever)
      as_tibble %>%
      select(sample(1:choose(nsurvive, nparents), nchild)) %>% # sample possible parent couples (as many as there will be children)
      as.list %>%
      lapply(
        function(ind) {
          parents[ind] %>%
            give_birth # spawn the children
        }
      )
    names(children) <- NULL
    
    # introduce some random variation
    mutants <- spawn_n(npop - nsurvive - nchild)
    
    # next round's population
    population <- c(parents, children, mutants)
  }
  
  # results --------------------------------------------------------------------
  cat("Completed", generation, "generations.")
  return(currentBest)
}
