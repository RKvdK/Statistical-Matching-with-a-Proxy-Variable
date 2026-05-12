# The SIM function generates a synthetic population 
# Including samples A, B, the overlap and external samples

SIM <- function(
    
  # ============================================================================
  # Function input parameters
  # ============================================================================
  
    nY = 3, # Number of Y categories
    nZ = 3, # Number of Z categories
    nP = 3, # Number of P categories
    nX = 6, # Number of X categories
    
    Ypopmar = rep(1/3, 3), # Vector of population margin Y
    Zpopmar = rep(1/3, 3), # Vector of population margin Z
    Xpopmar = rep(1/6, 6), # Vector of population margin X
    
    N = 1000000, # Population size
    
    tranmat_diag, # Strength relation Y and P 
    tranmat_sym = TRUE, # Symmetrical property of transition matrix
    
    tol = 1e-6, # Tolerance value
    seed = NULL, # Random seed for reproducibility
    
    n, # Sample size for samples A and B
    nE, # External sample size
    PercOverlap, # Sample size overlap percentage
    
    diagprob = 0.4, # Diagonal transition probability
    cia = 0, # If 0, CIA is valid
    
    mechanism = c("MAR", "MNAR"), # Selectivity mechanism
    pattern, # Selectivity pattern when mechanism is MAR
    degree # Selectivity degree when mechanism is MAR
    
) {
  
  # ============================================================================
  # Pre-requisite sanity checks
  # ============================================================================

  # Check whether population margins are positive
  
  if (any(Ypopmar < 0) || any(Zpopmar < 0) || any(Xpopmar < 0))
    stop("Population margins contain negative values.") 
  
  # Check whether population margins add up to 1
  
  if (abs(sum(Ypopmar) - 1) > tol) 
    stop("Sum of Ypopmar does not add up to 1.")
  if (abs(sum(Zpopmar) - 1) > tol)
    stop("Sum of Zpopmar does not add up to 1.")
  if (abs(sum(Xpopmar) - 1) > tol)
    stop("Sum of Xpopmar does not add up to 1.")
  
  # Check whether population margins align with the number of categories
  
  if (length(Ypopmar) != nY) 
    stop("Length of Ypopmar does not match nY.")
  if (length(Zpopmar) != nZ) 
    stop("Length of Zpopmar does not match nZ.")
  if (length(Xpopmar) != nX) 
    stop("Length of Xpopmar does not match nX.")
  
  # Check whether the number of categories across Y and P are equal
  
  if (nY != nP)
    stop("Number of categories across Y and P must be equal.")
  
  # ============================================================================
  # Input preparation
  # ============================================================================
  
  # A seed is only set when explicitly provided
  # Reproducibility across simulation conditions is handled in RUN()
   
  if (!is.null(seed)) set.seed(seed)
  
  # Match selectivity mechanism
  
  mechanism <- match.arg(as.character(mechanism), c("MAR", "MNAR")) 
  
  # Create labels for the variable levels
  
  Ylvl <- paste0("Y", seq_len(nY)) # Labels of Y levels
  Zlvl <- paste0("Z", seq_len(nZ)) # Labels of Z levels
  Xlvl <- paste0("X", seq_len(nX)) # Labels of X levels
  Plvl <- paste0("P", seq_len(nP)) # Labels of P levels
  
  # ============================================================================
  # Generate population Y variable
  # ============================================================================
  
  Ypop <- sample.int(nY, N, replace = TRUE, prob = Ypopmar)
  
  # ============================================================================
  # Generate population X variable
  # ============================================================================
  
  Xpop <- sample.int(nX, N, replace = TRUE, prob = Xpopmar)
  
  # ============================================================================
  # Generate population Y* variable
  # ============================================================================
  
  # Ensure that the diagonal of transition matrix same length as nY
  # If tranmat_diag is a single value, this value is replicated across nY
  
  if (length(tranmat_diag) == 1L) tranmat_diag <- rep(tranmat_diag, nY) 
  
  # Create memory for transition matrix
  
  tranmat <- matrix(0, nrow = nY, ncol = nP, dimnames = list(Y = Ylvl, P = Plvl))
  
  if (tranmat_sym) { # If transition matrix is symmetric
    
    off_diag <- (1 - tranmat_diag) / (nY - 1) # Define off-diagonal elements
    
    # Insert off-diagonal elements in transition matrix
    
    tranmat  <- matrix(rep(off_diag, each = nP), nrow = nY, byrow = TRUE)
    
    # Insert diagonal elements in transition matrix
    
    diag(tranmat) <- tranmat_diag
    dimnames(tranmat) <- list(Y = Ylvl, P = Plvl) # Ensure row and column names
    
  } else { # If transition matrix is not symmetric
    
    for (i in seq_len(nY)) {
      
      d <- tranmat_diag[i] # Define diagonal element
      massoff <- 1 - d # Off-diagonal probability mass                
      
      j <- i # Connect Y and P levels                    
      offcol <- setdiff(seq_len(nP), j) # Select all columns except the diagonal one
      
      if (length(offcol) == 0L) { # In case nP == 1
        
        tranmat[i, j] <- d # Only diagonal element
        
      } else { # In case nP > 1
        
        # Random off-diagonal elements scaled to off-diagonal mass
        
        off <- randoff(length(offcol)) * massoff
        
        tranmat[i, offcol] <- off   # Insert off-diagonal elements
        tranmat[i, j] <- d # Insert diagonal element
       
      }
    }
  }  
  
  # Memory vector for the proxy variable P
  
  Ppop <- integer(N) 
  
  for (i in seq_len(nY)) {
    
    idx <- (Ypop == i) # Indices for Y = i
    ni  <- sum(idx) # Number of units with Y = i
    
    # Generate P based on the transition matrix
    
    Ppop[idx] <- sample.int(nP, size = ni, replace = TRUE, prob = tranmat[i, ])
      
  }
  
  # ============================================================================
  # Generate population Z variable
  # ============================================================================
  
  Zpop <- genZ(Ypop, Ppop, Ylvl, Zlvl, Plvl, diagprob = diagprob, cia = cia) 
  
  # ============================================================================
  # Generate population ID variable
  # ============================================================================
  
  ID <- seq_len(N) 
  
  # ============================================================================
  # Sample A
  # ============================================================================
  
  idxA <- sample.int(N, n, replace = FALSE) # Sample A indices
  
  A <- data.frame(
    ID = idxA,
    Y  = Ypop[idxA],
    X  = Xpop[idxA]
  )
  
  # ============================================================================  
  # Sample B
  # ============================================================================
    
  Poverlap  <- PercOverlap # Overlap sample fractional size
  Noverlap  <- round(Poverlap * n) # Overlap sample size
  
  IDoverlap <- sample(idxA, Noverlap, replace = FALSE) # Overlap sample indices
  
  restB <- n - Noverlap # Remaining sample B size
  
  notA <- rep(TRUE, N); notA[idxA] <- FALSE # Population units not in sample A
  pool <- which(notA) # Pool of units outside sample A
  
  if (restB > length(pool)) stop("Not enough population units outside A for requested overlap.")
  
  subB <- sample(pool, restB, replace = FALSE) # Sample non-overlap indices
  idxB <- c(IDoverlap, subB) # Sample B indices
  
  B <- data.frame(
    ID = idxB,
    Z  = Zpop[idxB],
    X  = Xpop[idxB],
    P  = Ppop[idxB]
  )
  
  # ============================================================================
  # Overlap sample
  # ============================================================================
  
  ov <- IDoverlap
  
  overlap <- data.frame(
    ID = ov,
    Y  = Ypop[ov],
    X  = Xpop[ov],
    Z  = Zpop[ov],
    P  = Ppop[ov]
  )
  
  # ============================================================================
  # Sample E
  # ============================================================================
  
  if (mechanism  == "MAR") {
    
    # Retrieve pattern and degree selectivity weights for X
    
    wX <- mar_patdeg[[pattern]][[degree]] 
    
    # Fallback to original if there is no degree specified
    
    if (is.null(wX)) wX <- mar_patdeg[[pattern]][["original"]] 
    
    # Convey the weights to the population data
    
    wE <- wX[Xpop]
    
  } else {
    
    # Retrieve the YZ matrix for the specified pattern
    
    matYZ <- mnar_YZ[[pattern]] 
    
    # Convey the weights to the population data
    
    wE <- matYZ[cbind(Ypop, Zpop)]
    
    }
  
  # Ensure weights are non-negative
  
  wE <- pmax(wE, 0)
  
  # SRS as fallback in case all weights are zero
  
  if (all(wE == 0)) wE[] <- 1
  
  # Retrieve normalized probabilities
  
  probE <- wE / sum(wE)
  
  # Sample E based on indices drawn according to the weights
  
  idxE <- sample.int(N, nE, replace = FALSE, prob = probE)
  
  E <- data.frame(
    ID = idxE,
    Y  = Ypop[idxE],
    P  = Ppop[idxE]
  )
  
  # ============================================================================
  # Marginal distribution of X
  # ============================================================================
  
  marX <- marSM(B$X, levels = Xlvl, a = 0.5)
  
  # ============================================================================
  # True population joint distribution of Y and Z
  # ============================================================================
  
  YZpop <- prop.table(table(Ypop, as.integer(Zpop))) 
  dimnames(YZpop) <- list(Y = Ylvl, Z = Zlvl)
  
  # ============================================================================
  # Output list
  # ============================================================================
  
  SIMout <- list( # Output list
    meta = list(
      Xsource = "B",
      mechanism = mechanism,
      pattern = pattern,
      degree = degree), 
    poplvl = list( # List variable levels
      Ylvl = Ylvl, 
      Zlvl = Zlvl,
      Xlvl = Xlvl,
      Plvl = Plvl
    ),
    tranmat = tranmat, # Transition Matrix
    popjoint = list(
      YZpop = YZpop,
      YPpop = prop.table(table(Ypop, Ppop)),
      PconY = prop.table(table(Ypop, Ppop), margin = 1)),
    samples = list( # Samples A, B, E and the overlap
      A = A,
      B = B,
      E = E, 
      overlap = overlap),
    marX = marX
  )
  
  return(SIMout) 
  
}