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
  recordOutput = F # if performance is saved to file should only fitness be
  # recorded (F), or alsso the output of FUN (T)
) {
  # Performs a simple genetic optimization for hyper-parameter tuning
  
  if (nchild + nsurvive > npop - 1) stop("Reduce children or survivors")
  if (nchild > choose(nsurvive, nparents)) stop("Too many children requested")
  
  # helpers --------------------------------------------------------------------
  spawn_n <- function(n) {
    # helper - spawns n draws from individual_generaotr
    n %>%
      seq %>%
      lapply(
        function() {
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
    population <- checkpoint %>%
      read_csv %>%
      arrange(fitness) %>%
      slice(1:npop)
    population <- population["param_" %in% colnames(population)]
    colnames(population) <- str_replace(colnames(population), "param_", "")
    skipFirst <- T
  }
  
  # main loop ------------------------------------------------------------------
  for (generation in seq(generations)) {
    cat("------- Training Generation", generation, "---------\n")
    print(reduce(population, rbind))
    
    if (!skipFirst) {
      # calculate fitness
      fitness <- population %>%
        sapply(
          function(individual) {
            individual %>%
              FUN %>%
              lossWrapper %>%
              checkWriter(individual) %>%
              loss_message
          }
        )
      
      if (min(fitness) < tol) {
        cat("Converged to target. Best parameters are: \n")
        cat(population[which.min(fitness)])
        break()
      }
      
      # cutoff for survival
      maxfitness <- sort(fitness, decreasing = F)[nsurvive]
      
      # choose survivors      
      parents <- population[fitness <= minfitness]
      parents <- parents[1:nsurvive] # in case of tie
      
    } else {
      parents <- population[1:nsurvive] # if performance is loaded from 
      #checkpoint no need to reavaluate fitness
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
    
    # introduce some random variation
    mutants <- spawn_n(npop - nsurvive - nchild)
    
    # next round's population
    population <- c(parents, children, mutants)
    
  }
  
}