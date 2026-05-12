# The MC function performs a Monte Carlo simulation to evaluate the performance of the three different estimators 
# It resamples the input data multiple times, computes the estimators for each resample, and calculates bias and root mean square error (RMSE) 
# The SIM function is requisite for this function to work properly

MC <- function(SIMout, # SIM output
               MCnum, # Number of simulation sets
               seed = NULL,
               keep_iter = TRUE) { 
  
  # ============================================================================
  # SIM output
  # ============================================================================

  # Samples
  
  A0 <- SIMout$samples$A # Sample A 
  B0 <- SIMout$samples$B # Sample B
  E0 <- SIMout$samples$E # Sample E 
  Overlap0 <- SIMout$samples$overlap # Overlap sample 
  
  # Population variable levels
  
  Ylvl <- SIMout$poplvl$Ylvl # Y levels
  Zlvl <- SIMout$poplvl$Zlvl # Z levels
  Xlvl <- SIMout$poplvl$Xlvl # X levels
  Plvl <- SIMout$poplvl$Plvl # Proxy levels
  
  # True population joint distribution of Y and Z 
  
  YZpop <- SIMout$popjoint$YZpop 
  
  # Marginal distribution of X
  
  marXb <- SIMout$marX
  marXb <- marXb / sum(marXb) # Normalize the marginal distribution of X
  
  # ============================================================================
  # Memory arrays
  # ============================================================================
  
  # Memory arrays for the estimates of the three methods across Monte Carlo iterations
  
  I <- length(Ylvl) # Number of Y levels
  J <- length(Zlvl) # Number of Z levels
  
  # Doubly robust estimator
  
  theta_DRE <- array(NA_real_, dim = c(MCnum, I, J), dimnames = list(iter = seq_len(MCnum), Ylvl, Zlvl))
  
  # Iterative Proportional Fitting estimator
  
  theta_IPF <- array(NA_real_, dim = c(MCnum, I, J), dimnames = list(iter = seq_len(MCnum), Ylvl, Zlvl))
  
  # CIA proxy estimator
  
  theta_PROXY <- array(NA_real_, dim = c(MCnum, I, J), dimnames = list(iter = seq_len(MCnum), Ylvl, Zlvl))
  
  # Convert the true population joint distribution of Y and Z to a matrix and then to an array for easier calculations
  
  YZpop <- as.matrix(YZpop) # Convert into matrix
  YZpop_arr <- array(YZpop, dim = dim(theta_DRE)) # Convert into array
  
  # ============================================================================
  # Initialize for-loop
  # ============================================================================
  
  # A seed is only set when explicitly provided
  # Reproducibility across simulation conditions is handled in RUN()
  
  if (!is.null(seed)) set.seed(seed)
  
  for (b in seq_len(MCnum)) { # Monte Carlo loop
    
    # ==========================================================================
    # Bootstrap the samples
    # ==========================================================================
    
    # Generate indices for resampling with replacement for each sample 
    
    iA <- sample.int(nrow(A0), nrow(A0), replace = TRUE) # Indices for sample A
    iB <- sample.int(nrow(B0), nrow(B0), replace = TRUE) # Indices for sample B
    iE <- sample.int(nrow(E0), nrow(E0), replace = TRUE) # Indices for sample E
    
    # Resample the samples based on the generated indices
    
    Ab <- A0[iA, , drop = FALSE] # Resampled sample A
    Bb <- B0[iB, , drop = FALSE] # Resampled sample B
    Eb <- E0[iE, , drop = FALSE] # Resampled sample E
    
    if (nrow(Overlap0) == 0) { # If there is no unit overlap
      
      # Retain empty overlap sample
      
      Ob <- Overlap0
      
    } else { # If there is unit overlap
      
      iO <- sample.int(nrow(Overlap0), nrow(Overlap0), replace = TRUE) # Generate indices
      Ob <- Overlap0[iO, , drop = FALSE] # Resample the overlap sample based on the generated indices
      
    }
    
    # ==========================================================================
    # Compute pre-requisite conditional probabilities
    # ==========================================================================
    
    # Conditional probability of Y given P from sample E
    
    YconP_E <- conSM(Eb$Y, Eb$P, given = "col",
                     rlv = Ylvl, clv = Plvl, a = 0.5)
    
    # Conditional probability of P given X from sample B
    
    PconX_B <- conSM(Bb$P, Bb$X, given = "col",
                     rlv = Plvl, clv = Xlvl, a = 0.5)
    
    # List of conditional probabilities of Z given P, computed separately by X level
    
    ZconPX_B <- vector("list", length(Xlvl)) # Memory list
    names(ZconPX_B) <- Xlvl # Assign names to the list based on X levels
    
    for (xi in seq_along(Xlvl)) {
      
      idx <- Bb$X == xi
      ZconPX_B[[xi]] <- conSM(Bb$Z[idx], Bb$P[idx], given = "col",
                              rlv = Zlvl, clv = Plvl, a = 0.5)
      
    }
    
    # Apply alignment function to ensure that the P levels are consistent across the conditional probabilities
    
    aligned  <- alignP(YconP_E, PconX_B, ZconPX_B, Plvl) 
    
    YconP_E <- aligned$YconP # P(Y|PROXY) from sample E
    PconX_B  <- aligned$PconX # P(PROXY|X) from sample B
    ZconPX_B <- aligned$ZconPX # P(Z|PROXY,X) from sample B
    
    # ==========================================================================
    # CIA proxy estimator
    # ==========================================================================
    
    # Memory matrix
    
    YZ_CIA_full <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) 
    
    for (xi in seq_along(Xlvl)) {
      
      ZconPx <- ZconPX_B[[xi]] # Extract per level of X
      wPx <- as.numeric(PconX_B[, xi]) * as.numeric(marXb[xi]) # Weights
      YZ_CIA_full <- YZ_CIA_full + tcrossprod(YconP_E, sweep(ZconPx, 2, wPx, "*")) # Update YZ_CIA_full
      
    }
    
    YZ_CIA_full <- YZ_CIA_full / sum(YZ_CIA_full) # Normalize the resulting matrix
    
    jointPROXY <- orderYZ(YZ_CIA_full, Ylvl = Ylvl, Zlvl = Zlvl) # Store CIA proxy estimate
    
    # ==========================================================================
    # CIA(Y*) Overlap model + observed overlap
    # ==========================================================================
    
    # This part assumes that the overlap sample reflects the population structure correctly
    
    if (nrow(Ob) == 0) { # If there is no unit overlap
      
      YZ_CIA_overlap <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl),
                               dimnames = list(Ylvl, Zlvl)) # Zero matrix            
      
      YZ_overlap_obs <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl),
                               dimnames = list(Ylvl, Zlvl)) # Zero matrix       
      
    } else { # If there is unit overlap
      
      # Conditional probability of P given X (overlap) with smoothing
      
      PconX_O <- conSM(Ob$P, Ob$X, given = "col",
                       rlv = Plvl, clv = Xlvl, a = 0.5)
      
      # Conditional probability of Z given P and X (per X-level) with smoothing
      
      ZconPX_O <- vector("list", length(Xlvl)) # Memory list
      names(ZconPX_O) <- Xlvl # Assign names to the list based on X levels
      
      for (xi in seq_along(Xlvl)) {
        
        idx <- Ob$X == xi
        ZconPX_O[[xi]] <- conSM(Ob$Z[idx], Ob$P[idx], given = "col",
                                rlv = Zlvl, clv = Plvl, a = 0.5)
        
      }
      
      # Align P levels across Y|PROXY (from E), PROXY|X (overlap) and Z|PROXY,X (overlap)
      
      aligned_O <- alignP(YconP_E, PconX_O, ZconPX_O, Plvl)
      YconP_O   <- aligned_O$YconP
      PconX_O   <- aligned_O$PconX
      ZconPX_O  <- aligned_O$ZconPX
      
      YZ_CIA_overlap <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl),
                               dimnames = list(Ylvl, Zlvl)) # Memory matrix
      
      for (xi in seq_along(Xlvl)) {
        
        ZconPxO <- ZconPX_O[[xi]] # Extract per level of X
        wPxO <- as.numeric(PconX_O[, xi]) * as.numeric(marXb[xi]) # Weights
        YZ_CIA_overlap <- YZ_CIA_overlap + tcrossprod(YconP_O, sweep(ZconPxO, 2, wPxO, "*")) # Update YZ_CIA_overlap
        
      }
      
      YZ_CIA_overlap <- YZ_CIA_overlap / sum(YZ_CIA_overlap) # Normalize the resulting matrix
      
      # Extract the observed joint distribution of Y and Z in the overlap sample
      
      YZ_overlap_obs <- prop.table(table( # Joint distribution of Y and Z in overlap sample
        asFactorLevels(Ob$Y, Ylvl),
        asFactorLevels(Ob$Z, Zlvl)
      )) 
    }
    
    # ==========================================================================
    # Doubly Robust Estimator
    # ==========================================================================
    
    jointDRE <- (orderYZ(YZ_overlap_obs, Ylvl, Zlvl) - 
                   orderYZ(YZ_CIA_overlap, Ylvl, Zlvl)) + 
                     orderYZ(YZ_CIA_full, Ylvl, Zlvl) 
    
    jointDRE <- jointDRE / sum(jointDRE) # Normalize the resulting matrix
    
    # ==========================================================================
    # Iterative Proportional Fitting
    # ==========================================================================
    
    jointIPFx <- array(0, dim = c(length(Ylvl), length(Zlvl), length(Xlvl)),
                       dimnames = list(Ylvl, Zlvl, Xlvl)) # Memory array
    
    # Estimate target marginals for Y and Z conditional on X
    
    YconX_A  <- conSM(Ab$Y, Ab$X, given = "col", rlv = Ylvl, clv = Xlvl, a = 0.5)
    ZconX_Bb <- conSM(Bb$Z, Bb$X, given = "col", rlv = Zlvl, clv = Xlvl, a = 0.5)
    
    for (xi in seq_along(Xlvl)) { # Loop per level of X
      
      rowtar <- as.numeric(YconX_A[, xi, drop = TRUE]) # Row targets
      coltar <- as.numeric(ZconX_Bb[, xi, drop = TRUE]) # Column targets
      
      idxx <- Ob$X == xi # Select overlap units per level of X
      
      if (any(idxx)) { # Seed matrix per level of X in case there is unit overlap
        
        # Joint distribution in overlap per level of X
        
        tabO <- table( 
          asFactorLevels(Ob$Y[idxx], Ylvl),
          asFactorLevels(Ob$Z[idxx], Zlvl)
        )
        
        seedIPF <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory matrix
        seedIPF[rownames(tabO), colnames(tabO)] <- as.numeric(tabO) # Insert seed matrix values based on unit overlap
        
      } else { # If there are no overlapping units
        
        # Seed matrix is based on the outer product of the target mariginals
        
        seedIPF <- outer(rowtar, coltar) # Compute the outer product
        dimnames(seedIPF) <- list(Ylvl, Zlvl) # Assign the correct dimension names
        
      }
      
      jointIPFx[, , xi] <- ipfSafe(seedIPF, rowtar, coltar)$matrix # Perform IPF per level of X
      
    }
    
    jointIPF <- matrix(0, nrow = length(Ylvl), ncol = length(Zlvl), dimnames = list(Ylvl, Zlvl)) # Memory matrix
    
    for (xi in seq_along(Xlvl)) { # Loop per level of X
      
      # IPF estimate as weighted average of IPF estimates per level of X
      
      jointIPF <- jointIPF + jointIPFx[, , xi] * as.numeric(marXb[xi]) 
      
    }
    
    jointIPF <- jointIPF / sum(jointIPF) # Normalize the resulting matrix
    
    # Store the estimates of the three methods for the current Monte Carlo iteration
    
    theta_DRE[b, , ] <- jointDRE
    theta_IPF[b, , ] <- jointIPF
    theta_PROXY[b, , ] <- jointPROXY
    
  } # End of the Monte Carlo loop
  
  # ==========================================================================
  # Compute bias and variance summaries
  # ==========================================================================
  
  # Mean estimate per cell across Monte Carlo iterations
  
  mean_theta_DRE <- apply(theta_DRE, c(2,3), mean)
  mean_theta_IPF <- apply(theta_IPF, c(2,3), mean)
  mean_theta_PROXY <- apply(theta_PROXY, c(2,3), mean)
  
  # Bias per cell across Monte Carlo iterations
  
  bias_cell_DRE <- mean_theta_DRE - YZpop
  bias_cell_IPF <- mean_theta_IPF - YZpop
  bias_cell_PROXY <- mean_theta_PROXY - YZpop
  
  # Mean absolute bias per cell across Monte Carlo iterations
  
  absbias_cell_DRE <- apply(abs(theta_DRE - YZpop_arr), c(2,3), mean)
  absbias_cell_IPF <- apply(abs(theta_IPF - YZpop_arr), c(2,3), mean)
  absbias_cell_PROXY <- apply(abs(theta_PROXY - YZpop_arr), c(2,3), mean)
  
  # Variance per cell across Monte Carlo iterations
  
  var_cell_DRE <- apply(theta_DRE, c(2,3), var)
  var_cell_IPF <- apply(theta_IPF, c(2,3), var)
  var_cell_PROXY <- apply(theta_PROXY, c(2,3), var)
  
  # Standard deviation per cell across Monte Carlo iterations
  
  sd_cell_DRE <- sqrt(var_cell_DRE)
  sd_cell_IPF <- sqrt(var_cell_IPF)
  sd_cell_PROXY <- sqrt(var_cell_PROXY)
  
  
  # ============================================================================
  # Store the results
  # ============================================================================
  
  # Output list
  
  MCout <- list( # Output list
    
    MCnum = MCnum, # Number of Monte Carlo draws
    
    pop = list(
      YZpop = YZpop,
      Ylvl = Ylvl,
      Zlvl = Zlvl
    ),
    
    DRE = list( # Doubly robust estimator results
      theta_DRE = if (keep_iter) theta_DRE else NULL,
      mean_theta_DRE = mean_theta_DRE,
      bias_cell_DRE = bias_cell_DRE,
      absbias_cell_DRE = absbias_cell_DRE,
      var_cell_DRE = var_cell_DRE,
      sd_cell_DRE = sd_cell_DRE
    ),
    
    IPF = list( # Iterative Proportional Fitting results
      theta_IPF = if (keep_iter) theta_IPF else NULL,
      mean_theta_IPF = mean_theta_IPF,
      bias_cell_IPF = bias_cell_IPF,
      absbias_cell_IPF = absbias_cell_IPF,
      var_cell_IPF = var_cell_IPF,
      sd_cell_IPF = sd_cell_IPF
    ),
    
    PROXY = list( # CIA proxy results
      theta_PROXY = if (keep_iter) theta_PROXY else NULL,
      mean_theta_PROXY = mean_theta_PROXY,
      bias_cell_PROXY = bias_cell_PROXY,
      absbias_cell_PROXY = absbias_cell_PROXY,
      var_cell_PROXY = var_cell_PROXY,
      sd_cell_PROXY = sd_cell_PROXY
    )
  )
  
  return(MCout)
  
}