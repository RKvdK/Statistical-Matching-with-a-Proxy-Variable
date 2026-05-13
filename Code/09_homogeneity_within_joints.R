# This script evaluates the homogeneity of bias and standard deviation within the joint distribution

# ==============================================================================
# Extend cell_data with cell-level summaries for standard deviation
# ==============================================================================

data <- MC50 %>%
  
  select(-ID, -N) %>% 
  
  mutate( # Calculate cell-level summary statistics for bias, variance, standard deviation, and RMSE
    
    # Doubly robust estimator
    
    DRE_biasmean = map_dbl(DRE, ~ mean(.x$absbias_cell)),
    DRE_varmean  = map_dbl(DRE, ~ mean(.x$var_cell)),
    DRE_sdmean   = map_dbl(DRE, ~ mean(sqrt(.x$var_cell))),
    DRE_rmsemean = map_dbl(DRE, ~ mean(sqrt(.x$bias_cell_DRE^2 + .x$var_cell_DRE))),
    
    # Iterative proportional fitting
    
    IPF_biasmean = map_dbl(IPF, ~ mean(.x$absbias_cell)),
    IPF_varmean  = map_dbl(IPF, ~ mean(.x$var_cell)),
    IPF_sdmean   = map_dbl(IPF, ~ mean(sqrt(.x$var_cell))),
    IPF_rmsemean = map_dbl(IPF, ~ mean(sqrt(.x$bias_cell_IPF^2 + .x$var_cell_IPF))),
    
    # CIA proxy estimator
    
    PROXY_biasmean = map_dbl(PROXY, ~ mean(.x$absbias_cell)),
    PROXY_varmean  = map_dbl(PROXY, ~ mean(.x$var_cell)),
    PROXY_sdmean   = map_dbl(PROXY, ~ mean(sqrt(.x$var_cell))),
    PROXY_rmsemean = map_dbl(PROXY, ~ mean(sqrt(.x$bias_cell_PROXY^2 + .x$var_cell_PROXY))),
    
    # RMSE difference measures
    
    DIF_IPF_mean = IPF_rmsemean - DRE_rmsemean,
    DIF_PROXY_mean = PROXY_rmsemean - DRE_rmsemean 
    
  )

cell_data <- data %>%
  mutate(
    
    # Bias on the cell level

    DRE_biasmean = map_dbl(DRE, ~ mean(.x$absbias_cell)),
    DRE_biassd   = map_dbl(DRE, ~ sd(.x$absbias_cell)),
    
    IPF_biasmean = map_dbl(IPF, ~ mean(.x$absbias_cell)),
    IPF_biassd   = map_dbl(IPF, ~ sd(.x$absbias_cell)),
    
    PROXY_biasmean = map_dbl(PROXY, ~ mean(.x$absbias_cell)),
    PROXY_biassd   = map_dbl(PROXY, ~ sd(.x$absbias_cell)),
    
    # Standard deviation on the cell level
    
    DRE_sdmean = map_dbl(DRE, ~ mean(sqrt(.x$var_cell))),
    DRE_sdsd   = map_dbl(DRE, ~ sd(sqrt(.x$var_cell))),
    
    IPF_sdmean = map_dbl(IPF, ~ mean(sqrt(.x$var_cell))),
    IPF_sdsd   = map_dbl(IPF, ~ sd(sqrt(.x$var_cell))),
    
    PROXY_sdmean = map_dbl(PROXY, ~ mean(sqrt(.x$var_cell))),
    PROXY_sdsd   = map_dbl(PROXY, ~ sd(sqrt(.x$var_cell)))
  )
# ==============================================================================
# Cell data preparation 
# ==============================================================================

cell_summary_dat <- cell_data %>%
  
  # Add a row identifier for each parameter combination
  
  mutate(row_id = row_number()) %>%
  
  # Select the variables
  
  select(
    row_id, mechanism,
    DRE_biasmean, DRE_biassd, DRE_sdmean, DRE_sdsd,
    IPF_biasmean, IPF_biassd, IPF_sdmean, IPF_sdsd,
    PROXY_biasmean, PROXY_biassd, PROXY_sdmean, PROXY_sdsd
  ) %>%
  
  # Convert into long format
  
  pivot_longer(
    cols = -c(row_id, mechanism),
    names_to = c("estimator", "quantity", "stat"),
    names_pattern = "(DRE|IPF|PROXY)_(bias|sd)(mean|sd)",
    values_to = "value"
  ) %>%
  
  # Convert into wide format
  
  pivot_wider(
    id_cols = c(row_id, mechanism, estimator, quantity),
    names_from = stat,
    values_from = value
  ) %>%
  
  # Recode estimator variable and ensure factor levels for both the estimator and mechanism variables
  
  mutate(
    estimator = recode(estimator,
                       DRE = "Doubly robust",
                       IPF = "IPF",
                       PROXY = "CIA proxy"
    ),
    estimator = factor(estimator, levels = c("Doubly robust", "IPF", "CIA proxy")),
    quantity = factor(
      quantity,
      levels = c("bias", "sd"),
      labels = c("Bias", "Standard deviation")
    ),
    mechanism = factor(mechanism, levels = c("MAR", "MNAR"))
  )
# =============================================================================
# Cell data plot
# =============================================================================

cell_summary_plot <- ggplot(
  cell_summary_dat,
  aes(
    x = mean,
    y = sd,
    color = estimator
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.8
  ) +
  facet_grid(quantity ~ mechanism, scales = "free") +
  
  # Logarithmic scales 
  
  scale_x_log10() + 
  scale_y_log10() +
  
  scale_color_viridis_d(
    option = "C", # Colorblind-friendly palette
    begin = 0.1,
    end = 0.9
  ) +
  labs(
    x = "Mean over cells",
    y = "SD over cells",
    color = "Estimator"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    strip.background = element_blank(),
    strip.text = element_text(size = 20, face = "plain"),
    panel.grid.minor = element_blank()
  )

cell_summary_plot