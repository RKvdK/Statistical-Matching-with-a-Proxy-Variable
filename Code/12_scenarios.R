# This script specifies the scenarios that are subsets of the main grid
# These scenarios are studied in further detail

# ==============================================================================
# A: MAR Baseline
# ==============================================================================

A1 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.70,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

A2 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

# ==============================================================================
# B: MAR measurement error structure
# ==============================================================================

B3 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.70,
    tranmat_sym == FALSE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

B4 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.90,
    tranmat_sym == FALSE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

# ==============================================================================
# C: MAR Overlap extremes
# ==============================================================================

C5 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

C6 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0.30,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

# ==============================================================================
# D: MAR small data
# ==============================================================================

D7 <- simgrid %>%
  
  filter(
    n == 1000,
    nE == 200,
    PercOverlap == 0.10,
    tranmat_diag == 0.70,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

D8 <- simgrid %>%
  
  filter(
    n == 1000,
    nE == 200,
    PercOverlap == 0.10,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

# =============================================================================
# E: MAR heavy CIA violation
# ==============================================================================

E9 <- simgrid %>%
  
  filter(
    
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.70,
    tranmat_sym == TRUE,
    cia == 0.75,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

E10 <- simgrid %>%
  
  filter(
    
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.75,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

# =============================================================================
# F: MAR External extremes
# ==============================================================================

F11 <- simgrid %>%
  
  filter(
    
    n == 10000,
    nE == 2000,
    PercOverlap == 0.10,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

F12 <- simgrid %>%
  
  filter(
    
    n == 10000,
    nE == 8000,
    PercOverlap == 0.10,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MAR",
    pattern == "MAR_LinDec",
    degree == "original"
  )

# =============================================================================
# G: MNAR Baseline
# ==============================================================================

G13 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.70,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MNAR",
    pattern == "MNAR_ClasInc",
    degree == "none"
  )

G14 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MNAR",
    pattern == "MNAR_ClasInc",
    degree == "none"
  )

# =============================================================================
# H: MNAR interaction scenarios
# =============================================================================

H15 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MNAR",
    pattern == "MNAR_ClasInc_WeakInt",
    degree == "none"
  )

H16 <- simgrid %>%
  
  filter(
    n == 10000,
    nE == 5000,
    PercOverlap == 0.10,
    tranmat_diag == 0.90,
    tranmat_sym == TRUE,
    cia == 0.25,
    mechanism == "MNAR",
    pattern == "MNAR_ClasInc_StrongInt",
    degree == "none"
  )

scen <- list(
  
  # MAR baseline
  
  A1 = A1,
  A2 = A2,
  
  # MAR Measurement Error Structures
  
  B3 = B3,
  B4 = B4,
  
  # MAR overlap sample extremes
  
  C5 = C5,
  C6 = C6,
  
  # MAR small data (IPF)
  
  D7 = D7,
  D8 = D8,
  
  # MAR heavy CIA violation
  
  E9 = E9,
  E10 = E10,
  
  # MAR external sample extremes
  
  F11 = F11,
  F12 = F12,
  
  # MNAR baseline
  
  G13 = G13,
  G14 = G14,
  
  # MNAR interaction scenarios
  
  H15 = H15,
  H16 = H16
)