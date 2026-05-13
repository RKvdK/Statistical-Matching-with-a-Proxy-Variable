# Produce visualizations to enhance interpretation of the results

# ==============================================================================
# Data preparation
# ==============================================================================

# Filter MAR data

mardata <- data %>%
  filter(mechanism == "MAR")

# Filter MNAR data

mnardata <- data %>%
  filter(mechanism == "MNAR")

# ==============================================================================
# Missing at Random Sensitivity plot
# ==============================================================================

# Data preparation for MAR sensitivity plot

marpat_boxdat <- mardata %>%
  
  # Select the variables
  
  select( 
    pattern,
    DRE_biasmean, DRE_sdmean,
    IPF_biasmean, IPF_sdmean,
    PROXY_biasmean, PROXY_sdmean
  ) %>%
  
  # Convert into long format
  
  pivot_longer(
    cols = -pattern,
    names_to = c("estimator", "measure"),
    names_pattern = "(DRE|IPF|PROXY)_(biasmean|sdmean)",
    values_to = "value"
  ) %>%
  
  # Recode estimator names
  
  mutate(
    estimator = recode(estimator,
                       DRE = "Doubly robust",
                       IPF = "IPF",
                       PROXY = "CIA proxy"
    ),
    
    # Enforce factor structure for both estimator and pattern variables
    
    estimator = factor(estimator, levels = c("Doubly robust", "IPF", "CIA proxy")),
    measure = factor(
      measure,
      levels = c("biasmean", "sdmean"),
      labels = c("Mean absolute bias", "Mean standard deviation")
    ),
    pattern = recode(pattern,
                     MAR_ExtInc = "Extreme increase",
                     MAR_LinDec = "Linear decrease",
                     MAR_Step   = "Step function",
                     MAR_Ushape = "U-shaped"
    ),
    pattern = factor(
      pattern,
      levels = c("Linear decrease", "U-shaped", "Step function", "Extreme increase")
    )
  )

# MAR sensitivity plot

mar_pattern_plot <- ggplot(
  marpat_boxdat,
  aes(
    x = pattern,
    y = value
  )
) +
  geom_boxplot( # Boxplot
    outlier.alpha = 0.15, 
    width = 0.7,
    linewidth = 1.1,
    color = "black",
    fill = "white"
  ) +
  facet_grid( # Facet by measure and estimator
    measure ~ estimator,
    scales = "free_y",
    switch = "y"
  ) +
  scale_y_log10() + # Logarithmic scale for y-axis
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 26) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 24, face = "plain"), 
    strip.placement = "outside",
    strip.text.y.left = element_text(
      size = 24,
      angle = 90,
      vjust = 0.5
    ),
    axis.title = element_text(size = 25),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(
      size = 16,
      angle = 20,
      hjust = 1
    ),
    panel.grid.minor = element_blank()
  )

mar_pattern_plot

# ==============================================================================
# Missing Not At Random Sensitivity plot
# ==============================================================================

# Data preparation for MNAR sensitivity plot

mnarpat_boxdat <- mnardata %>%
  
  # Select the variables
  
  select(
    pattern,
    DRE_biasmean, DRE_sdmean,
    IPF_biasmean, IPF_sdmean,
    PROXY_biasmean, PROXY_sdmean
  ) %>%
  
  # Convert into long format
  
  pivot_longer(
    cols = -pattern,
    names_to = c("estimator", "measure"),
    names_pattern = "(DRE|IPF|PROXY)_(biasmean|sdmean)",
    values_to = "value"
  ) %>%
  
  # Recode estimator names
  
  mutate(
    estimator = recode(estimator,
                       DRE = "Doubly robust",
                       IPF = "IPF",
                       PROXY = "CIA proxy"
    ),
    
    # Enforce factor structure for both estimator and pattern variables
    
    estimator = factor(estimator, levels = c("Doubly robust", "IPF", "CIA proxy")),
    measure = factor(
      measure,
      levels = c("biasmean", "sdmean"),
      labels = c("Mean absolute bias", "Mean standard deviation")
    ),
    pattern = recode(pattern,
                     MNAR_ClasInc = "Classic increase",
                     MNAR_ClasInc_WeakInt = "Weak interaction",
                     MNAR_ClasInc_ModInt = "Moderate interaction",
                     MNAR_ClasInc_StrongInt = "Strong interaction",
                     MNAR_ClasInc_ExtInt = "Extreme interaction",
                     MNAR_NonMono = "Non-monotonic",
                     MNAR_Yonly = "Y only"
    ),
    pattern = factor(
      pattern,
      levels = c(
        "Classic increase",
        "Weak interaction",
        "Moderate interaction",
        "Strong interaction",
        "Extreme interaction",
        "Non-monotonic",
        "Y only"
      )
    )
  )

# MNAR sensitivity plot

mnar_pattern_plot <- ggplot(
  mnarpat_boxdat,
  aes(
    x = pattern,
    y = value
  )
) +
  geom_boxplot( # Boxplot
    outlier.alpha = 0.15,
    width = 0.7,
    linewidth = 1.1,
    color = "black",
    fill = "white"
  ) +
facet_grid( # Facet by measure and estimator
    measure ~ estimator, 
    scales = "free_y",
    switch = "y"
  ) +
  scale_y_log10() + # Logarithmic scale for y-axis
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 26) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 24, face = "plain"),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      size = 24,
      angle = 90,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(
      size = 15,
      angle = 20,
      hjust = 1
    ),
    
    panel.grid.minor = element_blank()
    
  )

mnar_pattern_plot

# ==============================================================================
# Driver plot data preparation
# ==============================================================================

driver_dat <- data %>%
  
  # Calculate the proportion of the external sample size relative to the total sample size
  
  mutate(
    nE_prop = nE / n) %>%  
  
  # Select the variables
  
  select(
    mechanism, n, nE_prop, PercOverlap, tranmat_diag, tranmat_sym, cia,
    DIF_IPF_mean, DIF_PROXY_mean
  ) %>%
  
  # Convert into long format
  
  pivot_longer(
    cols = c(DIF_IPF_mean, DIF_PROXY_mean),
    names_to = "comparison",
    values_to = "rmse_difference"
  ) %>%
  
  # Recode estimator names
  
  mutate(
    comparison = recode(
      comparison,
      DIF_IPF_mean = "IPF",
      DIF_PROXY_mean = "CIA proxy"
    ),
    
    # Enforce correct variable structures
    
    comparison = factor(comparison, levels = c("CIA proxy", "IPF")),
    n = factor(n),
    nE_prop = factor(nE_prop),
    PercOverlap = factor(PercOverlap),
    tranmat_diag = factor(tranmat_diag),
    tranmat_sym = factor(
      tranmat_sym,
      levels = c(TRUE, FALSE),
      labels = c("Symmetric", "Asymmetric")
    ),
    cia = factor(cia)
  ) %>%
  
  # Convert into long format
  
  pivot_longer(
    cols = c(n, nE_prop, PercOverlap, tranmat_diag, tranmat_sym, cia),
    names_to = "condition",
    values_to = "level"
  ) %>%
  
  # Recode condition names
  
  mutate(
    condition = recode(
      condition,
      n = "Sample size A/B",
      nE_prop = "External sample proportion",
      PercOverlap = "Unit overlap",
      tranmat_diag = "Proxy strength",
      tranmat_sym = "Error structure",
      cia = "CIA violation"
    ),
    
    # Enforce factor structure for condition variable
    
    condition = factor(
      condition,
      levels = c(
        "Sample size A/B",
        "External sample proportion",
        "Unit overlap",
        "Proxy strength",
        "Error structure",
        "CIA violation"
      )
    )
  )

# ==============================================================================
# MAR driver plot
# ==============================================================================

driver_plot_mar <- driver_dat %>%
  filter(mechanism == "MAR") %>%
  ggplot(aes(x = level, y = rmse_difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) + # Reference zero line
  geom_boxplot(
    outlier.alpha = 0.15,
    width = 0.65,
    linewidth = 1.1
  ) +
  facet_grid(
    comparison ~ condition, # Facet by comparison and condition
    
    # Use free scales and space for x-axis to allow for different numbers of levels across conditions
    # Preventing excessive spacing for conditions with fewer levels
    
    scales = "free_x",
    space = "free_x", #
    
    # Label modifications for better readability
    
    labeller = labeller(
      condition = c(
        "External sample proportion" = "External\nsample proportion",
        "Error structure" = "Error\nstructure"
      )
    )
  ) +
  labs(
    x = NULL,
    y = "RMSE difference (comparison estimator - DRE)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(
      size = 18,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 14,
      face = "plain",
      lineheight = 0.9
    ),
    panel.grid.minor = element_blank()
  )

driver_plot_mar

# ==============================================================================
# MNAR driver plot
# ==============================================================================

driver_plot_mnar <- driver_dat %>%
  filter(mechanism == "MNAR") %>%
  ggplot(aes(x = level, y = rmse_difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) + # Reference zero line
  geom_boxplot(
    outlier.alpha = 0.15,
    width = 0.65,
    linewidth = 1.1
  ) +
  facet_grid(
    comparison ~ condition, # Facet by comparison and condition
    
    # Use free scales and space for x-axis to allow for different numbers of levels across conditions
    # Preventing excessive spacing for conditions with fewer levels
    
    scales = "free_x",
    space = "free_x",
    
    # Label modifications for better readability
    
    labeller = labeller(
      condition = c(
        "External sample proportion" = "External\nsample proportion",
        "Error structure" = "Error\nstructure"
      )
    )
  ) +
  labs(
    x = NULL,
    y = "RMSE difference (comparison estimator - DRE)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(
      size = 18,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 16),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 14,
      face = "plain",
      lineheight = 0.9
    ),
    panel.grid.minor = element_blank()
  )

driver_plot_mnar

# ==============================================================================
# Step function data preparation
# ==============================================================================

step_dat <- data %>%
  
  # Calculate absolute external sample size and absolute unit overlap
  
  mutate(
    n_group = factor(n, levels = sort(unique(n))),
    nE_abs = nE,
    nAB_abs = n * PercOverlap
    
  ) %>%
  
  # Select the variables
  
  select(
    mechanism, n_group, nE_abs, nAB_abs,
    DIF_IPF_mean, DIF_PROXY_mean
  ) %>%
  
  # Convert into long format 
  
  pivot_longer(
    cols = c(DIF_IPF_mean, DIF_PROXY_mean),
    names_to = "comparison",
    values_to = "rmse_difference"
  ) %>%
  
  # Recode variable names
  
  mutate(
    comparison = recode(
      comparison,
      DIF_IPF_mean = "IPF",
      DIF_PROXY_mean = "CIA proxy"
    ),
    
    # Enforce factor structure for comparison variable
    
    comparison = factor(comparison, levels = c("CIA proxy", "IPF"))
  ) %>%
  
  # Convert into long format
  
  pivot_longer(
    cols = c(nE_abs, nAB_abs),
    names_to = "condition",
    values_to = "absolute_size"
  ) %>%
  
  # Recode condition variable
  mutate(
    condition = recode(
      condition,
      nE_abs = "External sample size",
      nAB_abs = "Unit overlap"
    ),
    
    # Enforce factor structure for condition and absolute size variables
    
    condition = factor(
      condition,
      levels = c("External sample size", "Unit overlap")
    ),
    absolute_size = factor(
      absolute_size,
      levels = sort(unique(absolute_size))
    )
  )

# ==============================================================================
# MAR stepfunction plot
# ==============================================================================

step_plot_mar <- step_dat %>%
  filter(mechanism == "MAR") %>%
  ggplot(aes(
    x = absolute_size,
    y = rmse_difference,
    fill = n_group,
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) + # Reference zero line
  geom_boxplot(
    outlier.alpha = 0.15,
    width = 0.65,
    linewidth = 1.1,
    position = position_dodge2(preserve = "single")
  ) + 
  
  # Facet by comparison and condition with free x-axis
  
  facet_grid(
    comparison ~ condition,
    
    # Free scales and space for x-axis to allow for different numbers of levels across conditions
    
    scales = "free_x",
    space = "free_x"
  ) + 
  
  scale_fill_manual( # Greyscale fills for sample size groups
    values = c(
      "1000" = "white",
      "10000" = "grey70",
      "100000" = "gray35"
    )
  ) +
  labs(
    x = NULL,
    y = "RMSE difference (comparison estimator - DRE)",
    fill = "Sample size A/B"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.7, "cm"),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "plain"),
    panel.grid.minor = element_blank()
  )

step_plot_mar

# ==============================================================================
# MNAR stepfunction plot
# ==============================================================================

step_plot_mnar <- step_dat %>%
  filter(mechanism == "MNAR") %>%
  ggplot(aes(x = absolute_size, y = rmse_difference, fill = n_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) + # Reference zero line
  geom_boxplot(
    outlier.alpha = 0.15,
    width = 0.65,
    linewidth = 1.1,
    position = position_dodge2(preserve = "single")
  ) +
  
  # Facet by comparison and condition with free x-axis
  
  facet_grid(
    comparison ~ condition,
    
    # Free scales and space for x-axis to allow for different numbers of levels across conditions
    
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual( # Greyscale fills for sample size groups
    values = c(
      "1000" = "white",
      "10000" = "grey70",
      "100000" = "gray35"
    )
  ) +
  labs(
    x = NULL,
    y = "RMSE difference (comparison estimator - DRE)",
    fill = "Sample size A/B"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.7, "cm"),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "plain"),
    panel.grid.minor = element_blank()
  )

step_plot_mnar