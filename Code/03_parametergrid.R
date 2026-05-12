# The GRID function creates a parameter grid for the Monte Carlo Simulation

# ==============================================================================
# Default grid
# ==============================================================================

GRID <- function(
    n = c(1000, 10000, 100000), # Size of samples A and B
    percExternal = c(0.20, 0.50, 0.80), # External sample size proportional to samples A and B
    percOverlap = c(0.00, 0.10, 0.30), # Overlap sample size proportional to samples A and B
    p_diag = c(0.50, 0.70, 0.90), # Association strength between Y and its proxy
    w  = c("Symmetrical","Asymmetrical"), # Whether the misclassification probabilities are balanced
    CIA = c(0, 0.25, 0.75), # Degree of violation of the conditional independence assumption
    N = 1e6, # Population size
    diagprob = 0.4, # Diagonal probability of the transition matrix
    id_start = 1L, # ID variable starting value
    nE_fun = function(n, percExt) round(n * percExt) # Function to compute nE
) {
  
  # ============================================================================
  # Missing at Random
  # ============================================================================
  
  # MAR patterns
  
  MAR_patterns <- c("MAR_LinDec", # Linear Decrease
                    "MAR_Ushape", # U-Shaped
                    "MAR_Step",  # Step Function
                    "MAR_ExtInc") # Extreme Increase
  
  
  # Mar degrees
  
  MAR_degrees <- c("reduced", # Reduced selectivity degree
                   "original", # Original selectivity degree
                   "increased") # Increased selectivity degree
  
  # MAR grid
  
  base_mar <- expand.grid(
    n = n,
    nEn = percExternal,
    nABn = percOverlap,
    p = p_diag,
    w = w,
    CIA = CIA,
    mechanism = "MAR",
    pattern = MAR_patterns,
    stringsAsFactors = FALSE
  ) %>%
    tidyr::crossing(degree = MAR_degrees) %>%
    filter(
      pattern == "MAR_LinDec" | degree == "original"
    )
  
  # ============================================================================
  # Missing Not at Random
  # ============================================================================
  
  # MNAR main and interaction scenarios
  
  MNAR_patterns <- c(
    
    # Main scenarios
    
    "MNAR_ClasInc",
    "MNAR_NonMono",
    "MNAR_Yonly",
    
    # Interaction scenarios
    
    "MNAR_ClasInc_WeakInt",
    "MNAR_ClasInc_ModInt",
    "MNAR_ClasInc_StrongInt",
    "MNAR_ClasInc_ExtInt"
  )
  
  # MNAR grid
  
  base_mnar <- expand.grid(
    n = n,
    nEn = percExternal,
    nABn  = percOverlap,
    p = p_diag,
    w = w,
    CIA  = CIA,
    mechanism = "MNAR",
    pattern = MNAR_patterns,
    degree = "none",
    stringsAsFactors = FALSE
  )
  
  # ============================================================================
  # Output grid
  # ============================================================================
  
  # Combine both grids
  
  base <- rbind(base_mar, base_mnar) 
  
  # Define grid output
  
  out <- transform(
    base,
    ID = seq.int(id_start, length.out = nrow(base)),
    tranmat_sym = (w == "Symmetrical"),
    tranmat_diag = p,
    PercOverlap = nABn,
    n = n,
    nE = nE_fun(n, nEn),
    N = N,
    cia = CIA,
    diagprob = diagprob
  )
  
  out <- out[, c("ID","N","n","nE","PercOverlap","tranmat_diag","tranmat_sym","cia","diagprob","mechanism","pattern","degree")] # Ensure variable order
  rownames(out) <- NULL # Ensure row names are clean
  
  out
  
}

simgrid <- GRID()
