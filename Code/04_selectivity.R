# This script contains all selectivity patterns and degrees assuming Missing at Random
# Also, this script contains all main and interaction scenarios assuming Missing Not at Random
# The derived patterns and scenarios are retrieved from Sojka (2025)

# ==============================================================================
# Missing at Random
# ==============================================================================

# Selectivity patterns and degrees

mar_patdeg <- list(
  MAR_LinDec = list( # Linear decrease pattern
    reduced = c(X1 = 3, X2 = 2.6, X3 = 2.2, X4 = 1.8, X5 = 1.4, X6 = 1), # Reduced selectivity degree
    original = c(X1 = 6, X2 = 5, X3 = 4, X4 = 3, X5= 2, X6 = 1), # Original selectivity degree
    increased = c(X1 = 12, X2 = 9.8, X3 = 7.6, X4 = 5.4, X5 = 3.2, X6 = 1) # Increased selectivity degree
    ),
  MAR_Ushape = list( # U-shaped pattern
    original = c(X1 = 8, X2 = 2, X3 = 1, X4 = 1, X5 = 2, X6 = 8) # Original selectivity degree
  ),
  MAR_Step = list( # Step function pattern
    original = c(X1 = 4, X2 = 4, X3 = 4, X4 = 1, X5 = 1, X6 = 1) # Original selectivity degree
  ),
  MAR_ExtInc = list( # Extreme increase pattern
    original = c(X1 = 1, X2 = 2, X3 = 3, X4 = 5, X5 = 10, X6 = 40) # Original selectivity degree
  )
)

# ==============================================================================
# MAR Pattern and degrees plots
# ==============================================================================

patdat <- imap_dfr(mar_patdeg, function(pattern_list, pattern_name) { # Iterate over patterns
  imap_dfr(pattern_list, function(values, degree_name) {# Iterate over degrees within patterns
    data.frame(
      pattern_raw = pattern_name, # Extract pattern name
      degree = degree_name, # Extract degree name
      X = as.numeric(gsub("X", "", names(values))),
      score = as.numeric(values) # Extract scores and convert to numeric
    )
  })
}) %>%
  mutate(
    pattern = recode(pattern_raw, # Recode raw pattern names
                     MAR_LinDec = "Linear Decrease",
                     MAR_Ushape = "U–Shaped",
                     MAR_Step   = "Step Function",
                     MAR_ExtInc = "Extreme Increase"
    ),
    pattern = factor( # Convert pattern to factor variable
      pattern,
      levels = c("Linear Decrease", "U–Shaped", "Step Function", "Extreme Increase")
    ),
    degree = factor( # Convert degree to factor variable
      degree,
      levels = c("reduced", "original", "increased")
    )
  ) 

# Pattern and degree plot

marpat <- ggplot(patdat, aes(x = X, y = score, group = degree)) +
  
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.7) + # Add reference line at score 1
  geom_point(aes(shape = degree), size = 3.5, color = "black") + # Add points for each score
  geom_line(aes(linetype = degree), linewidth = 1.4, color = "black") + # Add lines connecting points
  facet_wrap(~ pattern, scales = "free_y", nrow = 1) + # Facet by pattern with free y-axis
  scale_x_continuous(breaks = 1:6) +  # Set x-axis breaks at 1 to 6
  scale_y_continuous( # Set y-axis breaks
    breaks = c(1, 2, 4, 6, 8, 10, 12, 20, 30, 40)
  ) +
  
  # Line types for degrees
  
  scale_linetype_manual( 
    values = c(
      reduced = "dotted",
      original = "solid",
      increased = "dashed"
    ),
    breaks = c("reduced", "original", "increased"), # Degree breaks
    labels = c("Reduced", "Original", "Increased") # Degree labels
  ) +
  
  # Point shapes for degrees
  
  scale_shape_manual( 
    values = c(
      reduced = 16, # Circles
      original = 17, # Triangles
      increased = 15 # Squares
    ),
    breaks = c("reduced", "original", "increased"), # Degree breaks
    labels = c("Reduced", "Original", "Increased") # Degree labels
  ) +
  
  labs( 
    x = expression(italic(X)), # Italicize x-axis label
    y = "Selection Score",
    linetype = NULL,
    shape = NULL
  ) +
  
  theme_minimal() + 
  theme(
    strip.text = element_text(size = 24), # Facet label size
    axis.title = element_text(size = 25), # Axis title size
    axis.text = element_text(size = 21), # Axis text size
    legend.text = element_text(size = 19), # Legend text size
    panel.grid.minor = element_line(color = "grey85"), # Minor grid lines
    panel.grid.major = element_line(color = "grey80"), # Major grid lines
    legend.position = "bottom" # Legend position at the bottom
  )

marpat

# ==============================================================================
# Missing Not at Random
# ==============================================================================

mnar_YZ <- list(
  
  # Main scenarios
  
  MNAR_ClasInc = matrix( # Classic increase main scenario
    c(
      1, 2, 4,
      2, 4, 8,
      3, 6, 12
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  MNAR_NonMono = matrix( # Non-monotonic main scenario
    c(
      6, 2, 12,
      15, 5, 30,
      3, 1, 6
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  MNAR_Yonly = matrix( # Y-only main scenario
    c(
      4, 4, 4,
      1, 1, 1,
      7, 7, 7
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  # Interaction scenarios
  
  MNAR_ClasInc_WeakInt = matrix( # Classic increase with weak interaction
    c(
      1, 1.7, 3.3,  
      1.7, 4.3, 6.7,  
      2.5, 5, 11  
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  MNAR_ClasInc_ModInt = matrix( # Classic increase with moderate interaction
    c(
      1, 1, 2,
      1, 5, 4,
      1.5, 3, 10.8
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),
  
  MNAR_ClasInc_StrongInt = matrix( # Classic increase with strong interaction
    c(
      1.9, 1, 1.7,
      1, 8.8, 4,
      1.3, 3, 18.8
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  ),

  MNAR_ClasInc_ExtInt = matrix( # Classic increase with extreme interaction
    c(
      10, 1, 1,
      1, 10, 1,
      1, 1, 10
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      Y = c("Y1","Y2","Y3"),
      Z = c("Z1","Z2","Z3")
    )
  )
)

