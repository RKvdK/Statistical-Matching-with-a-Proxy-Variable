# The RUN function executes the simulation study over a grid of parameter settings
# The SIM and MC functions are requisite for this function to work

RUN <- function(fungrid, MCnum, keep_iter = FALSE) { 
  
  # ============================================================================
  # Initialize parallelisation
  # ============================================================================
  
  t0 <- Sys.time() # Start time
  
  n_workers <- max(1, parallel::detectCores() - 1) # Define the number of workers
  
  cl <- parallel::makeCluster(n_workers, outfile = 'MC50.txt') # Create cluster
  
  on.exit(parallel::stopCluster(cl), add = TRUE) # Stop cluster when RUN exits 
  
  parallel::clusterSetRNGStream(cl, iseed = 123) # Reproducible Random Number Generating 
  
  parallel::clusterExport(cl, # Load functions for workers
                          varlist = c("SIM", "MC", "conSM", "marSM",
                                      "alignP", "orderYZ", "asFactorLevels",
                                      "ipfSafe", "genZ", "randoff", "mar_patdeg", "mnar_YZ"),
                          envir = environment())
  
  outlist <- parallel::parLapply(cl, # Run the simulation in parallel over the parameter grid
                                 X = seq_len(nrow(fungrid)),
                                 fun = function(i){
    
    message(Sys.time(), ": started condition ", i)
    
    library(tidyverse)
    library(ggplot2)
                                   
  # ============================================================================
  # Select conditions rowwise from parameter grid
  # ============================================================================
    
    g <- fungrid[i, , drop = FALSE]  # Select row of the parameter grid
    
    S <- SIM( 
      N = g$N, # Population size
      n = g$n, # Sample size of samples A and B   
      nE = g$nE, # External sample size      
      PercOverlap = g$PercOverlap, # Percentage overlap between samples A and B
      tranmat_diag = g$tranmat_diag, # Diagonal elements of the transition matrix
      tranmat_sym = g$tranmat_sym, # Whether the transition matrix is symmetric 
      seed = as.integer(g$ID[[1]]), # Random seed for reproducibility       
      diagprob = g$diagprob, # Diagonal probability of the transition matrix      
      cia = g$cia, # Conditional independence assumption parameter    
      
      mechanism = if ("mechanism" %in% names(g)) as.character(g$mechanism[[1]]) else "MAR", # Selectivity mechanism
      pattern = if ("pattern" %in% names(g)) as.character(g$pattern[[1]]) else NA_character_, # Selectivity pattern
      degree = if ("degree" %in% names(g)) as.character(g$degree[[1]]) else NA_character_ # Selectivity degree
    )
    
    # ==========================================================================
    # Perform Monte Carlo
    # ==========================================================================

    MCsim <- MC(S, MCnum = MCnum, seed = g$ID, keep_iter = keep_iter) # Monte Carlo simulation 
    
    # ==========================================================================
    # Store results
    # ==========================================================================
    
    return(data.frame(
      
      ID = as.integer(g$ID), # Simulation ID
      N = as.integer(g$N), # Population size
      n = as.integer(g$n), # Sample size of samples A and B
      nE = as.integer(g$nE), # External sample size
      PercOverlap = as.numeric(g$PercOverlap), # Percentage overlap between samples A and B
      tranmat_diag = as.numeric(g$tranmat_diag), # Diagonal elements of the transition matrix
      tranmat_sym = g$tranmat_sym, # Whether the transition matrix is symmetric
      cia = as.numeric(g$cia), # Conditional independence assumption parameter
      diagprob = as.numeric(g$diagprob), # Diagonal probability of the transition matrix
      
      mechanism = if ("mechanism" %in% names(g)) as.character(g$mechanism) else "MAR", # Selectivity mechanism
      pattern = if ("pattern" %in% names(g)) as.character(g$pattern) else NA_character_, # Selectivity pattern
      degree = if ("degree" %in% names(g)) as.character(g$degree) else NA_character_, # Selectivity degree
      
      # Estimator results
      
      DRE = I(list(MCsim$DRE)),
      IPF = I(list(MCsim$IPF)),
      PROXY = I(list(MCsim$PROXY))
      
    ))
  })
  
  t1 <- Sys.time()  # End time
  message("RUN finished in: ", format(t1 - t0))
  
  return(do.call(rbind, outlist)) # Combine all list elements
  
}