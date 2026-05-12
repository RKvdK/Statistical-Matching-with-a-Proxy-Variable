# This script performs the simulations for the specified situations

# ==============================================================================
# Scenarios simulation
# ==============================================================================

scengrid <- rbind(A1, A2, B3, B4, C5, C6, D7, D8, E9, E10, F11, F12, G13, G14, H15, H16)

# scenres <- RUN(scengrid, MCnum = 5000)

# saveRDS(scenres, "Scenario Calculation Results")

scenres <- readRDS("Scenario Calculation Results")

# view(scenres)

# ==============================================================================
# Results visualization
# ==============================================================================

scenario_names <- names(scen) # Create scenario name vector

final_estimates <- data.frame(
  
  # Scenario names
  
  scenario = scenario_names,
  
  # Scenario-level bias values
  
  bias_DRE = sapply(scenres$DRE, function(x) mean(x$bias_cell_DRE, na.rm = TRUE)),
  bias_IPF = sapply(scenres$IPF, function(x) mean(x$bias_cell_IPF, na.rm = TRUE)),
  bias_PROXY = sapply(scenres$PROXY, function(x) mean(x$bias_cell_PROXY, na.rm = TRUE)),
  
  # Scenario-level absolute bias values
  
  absbias_DRE = sapply(scenres$DRE, function(x) mean(abs(x$bias_cell_DRE), na.rm = TRUE)),
  absbias_IPF = sapply(scenres$IPF, function(x) mean(abs(x$bias_cell_IPF), na.rm = TRUE)),
  absbias_PROXY = sapply(scenres$PROXY, function(x) mean(abs(x$bias_cell_PROXY), na.rm = TRUE)),
  
  # Scenario-level standard deviation values
  
  sd_DRE = sapply(scenres$DRE, function(x) mean(x$sd_cell_DRE, na.rm = TRUE)),
  sd_IPF = sapply(scenres$IPF, function(x) mean(x$sd_cell_IPF, na.rm = TRUE)),
  sd_PROXY = sapply(scenres$PROXY, function(x) mean(x$sd_cell_PROXY, na.rm = TRUE))
)

# Scenario names as factors with the same order as in the data frame

final_estimates$scenario <- factor(final_estimates$scenario, levels = final_estimates$scenario)

# Scenario-level MSE values

final_estimates$mse_DRE <- sapply(scenres$DRE, function(x) {
  mean(x$bias_cell_DRE^2 + x$var_cell_DRE, na.rm = TRUE)
})
  
final_estimates$mse_IPF <- sapply(scenres$IPF, function(x) {
  mean(x$bias_cell_IPF^2 + x$var_cell_IPF, na.rm = TRUE)
})

final_estimates$mse_PROXY <- sapply(scenres$PROXY, function(x) {
  mean(x$bias_cell_PROXY^2 + x$var_cell_PROXY, na.rm = TRUE)
})

# Scenario-level RMSE values

final_estimates$rmse_DRE <- sqrt(final_estimates$mse_DRE)
final_estimates$rmse_IPF <- sqrt(final_estimates$mse_IPF)
final_estimates$rmse_PROXY <- sqrt(final_estimates$mse_PROXY)


# saveRDS(final_estimates, "Final Estimates")

# ==============================================================================
# MAR plot
# ==============================================================================

# Filter MAR scenario

marest <- final_estimates %>%
  filter(!scenario %in% c("G13", "G14", "H15", "H16"))

# MAR absolute bias plot

marabsbias <- barplotSM(
  marest,
  cols = c("absbias_DRE", "absbias_IPF", "absbias_PROXY"),
  value_name = "Absolute Bias",
  title = "Absolute bias per scenario"
)

mar_absbias <- marabsbias +
  scale_fill_manual(
    name = "Estimator type",
    values = c(
      "DRE" = "white",
      "IPF" = "grey70",
      "PROXY" = "grey35"
    ),
    labels = c(
      "DRE" = "Doubly Robust",
      "IPF" = "Iterative Proportional Fitting",
      "PROXY" = "CIA proxy estimator"
    )
  )

# MAR standard deviation plot

marsd <- barplotSM(
  marest,
  cols = c("sd_DRE", "sd_IPF", "sd_PROXY"),
  value_name = "Standard Deviation",
  title = "Standard deviation per scenario"
)

mar_sd <- marsd +
  scale_fill_manual(
    name = "Estimator type",
    values = c(
      "DRE" = "white",
      "IPF" = "grey70",
      "PROXY" = "grey35"
    ),
    labels = c(
      "DRE" = "Doubly Robust",
      "IPF" = "Iterative Proportional Fitting",
      "PROXY" = "CIA proxy estimator"
    )
  )

# MAR RMSE plot

marrmse <- barplotSM(
  marest,
  cols = c("rmse_DRE", "rmse_IPF", "rmse_PROXY"),
  value_name = "RMSE",
  title = "RMSE per scenario"
)

mar_rmse <- marrmse +
  scale_fill_manual(
    name = "Estimator type",
    values = c(
      "DRE" = "white",
      "IPF" = "grey70",
      "PROXY" = "grey35"
    ),
    labels = c(
      "DRE" = "Doubly Robust",
      "IPF" = "Iterative Proportional Fitting",
      "PROXY" = "CIA proxy estimator"
    )
  )

# MAR excluding small samples sizes context

marest_filt <- final_estimates %>%
  filter(!scenario %in% c("D7", "D8", "G13", "G14", "H15", "H16"))

# MAR filtered absolute bias plot

marabsbias_filt <- barplotSM(
  marest_filt,
  cols = c("absbias_DRE", "absbias_IPF", "absbias_PROXY"),
  value_name = "Absolute Bias",
  title = "Absolute bias per scenario"
)

mar_absbias_filt <- marabsbias_filt +
  scale_fill_manual(
    name = "Estimator type",
    values = c(
      "DRE" = "white",
      "IPF" = "grey70",
      "PROXY" = "grey35"
    ),
    labels = c(
      "DRE" = "Doubly Robust",
      "IPF" = "Iterative Proportional Fitting",
      "PROXY" = "CIA proxy estimator"
    )
  )

# MAR filtered standard deviation plot

marsd_filt <- barplotSM(
  marest_filt,
  cols = c("sd_DRE", "sd_IPF", "sd_PROXY"),
  value_name = "Standard Deviation",
  title = "Standard deviation per scenario"
)

mar_sd_filt <- marsd_filt +
  scale_fill_manual(
    name = "Estimator type",
    values = c(
      "DRE" = "white",
      "IPF" = "grey70",
      "PROXY" = "grey35"
    ),
    labels = c(
      "DRE" = "Doubly Robust",
      "IPF" = "Iterative Proportional Fitting",
      "PROXY" = "CIA proxy estimator"
    )
  )

# MAR filtered RMSE plot

marrmse_filt <- barplotSM(
  marest_filt,
  cols = c("rmse_DRE", "rmse_IPF", "rmse_PROXY"),
  value_name = "RMSE",
  title = "RMSE per scenario"
)

mar_rmse_filt <- marrmse_filt +
  scale_fill_manual(
    name = "Estimator type",
    values = c(
      "DRE" = "white",
      "IPF" = "grey70",
      "PROXY" = "grey35"
    ),
    labels = c(
      "DRE" = "Doubly Robust",
      "IPF" = "Iterative Proportional Fitting",
      "PROXY" = "CIA proxy estimator"
    )
  )

# Construct MAR integrated figure

# Unfiltered MAR plots

mar_absbias_left <- mar_absbias +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18)
  )

mar_sd_left <- mar_sd +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18)
  )

mar_rmse_left <- mar_rmse +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

# Filtered MAR plots

mar_absbias_right <- mar_absbias_filt +
  scale_y_continuous(position = "right") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

mar_sd_right <- mar_sd_filt +
  scale_y_continuous(position = "right") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

mar_rmse_right <- mar_rmse_filt +
  scale_y_continuous(position = "right") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18)
  )

# Construct comprehensive MAR figure

marplot_clean <- (
  (mar_absbias_left | mar_absbias_right) /
    (mar_sd_left     | mar_sd_right) /
    (mar_rmse_left   | mar_rmse_right)
) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.box.margin = margin(t = 6)
  )

marplot_clean

# ==============================================================================
# MNAR plot
# ==============================================================================

# Filter MNAR scenarios

mnarest <- final_estimates %>%
  filter(scenario %in% c("G13", "G14", "H15", "H16"))

# MNAR absolute bias plot

mnar_absbias <- barplotSM(
  mnarest,
  cols = c("absbias_DRE", "absbias_IPF", "absbias_PROXY"),
  value_name = "Absolute Bias",
  title = "Absolute bias per scenario"
) +
  scale_fill_manual(
    name = "Estimator type",
    values = c(
      "DRE" = "white",
      "IPF" = "grey70",
      "PROXY" = "grey35"
    ),
    labels = c(
      "DRE" = "Doubly Robust",
      "IPF" = "Iterative Proportional Fitting",
      "PROXY" = "CIA proxy estimator"
    )
  )

mnar_absbias_clean <- mnar_absbias +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 15)
  )

# MNAR standard deviation plot

mnar_sd <- barplotSM(
  mnarest,
  cols = c("sd_DRE", "sd_IPF", "sd_PROXY"),
  value_name = "Standard Deviation",
  title = "Standard deviation per scenario"
) +
  scale_fill_manual(
    name = "Estimator type",
    values = c(
      "DRE" = "white",
      "IPF" = "grey70",
      "PROXY" = "grey35"
    ),
    labels = c(
      "DRE" = "Doubly Robust",
      "IPF" = "Iterative Proportional Fitting",
      "PROXY" = "CIA proxy estimator"
    )
  )

mnar_sd_clean <- mnar_sd +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 15)
  )

# Filter MNAR RMSE plot

mnar_rmse <- barplotSM(
  mnarest,
  cols = c("rmse_DRE", "rmse_IPF", "rmse_PROXY"),
  value_name = "RMSE",
  title = "RMSE per scenario"
) +
  scale_fill_manual(
    name = "Estimator type",
    values = c(
      "DRE" = "white",
      "IPF" = "grey70",
      "PROXY" = "grey35"
    ),
    labels = c(
      "DRE" = "Doubly Robust",
      "IPF" = "Iterative Proportional Fitting",
      "PROXY" = "CIA proxy estimator"
    )
  )

mnar_rmse_clean <- mnar_rmse +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )

# Construct comprehensive MNAR figure

mnarplot_clean <- (
  mnar_absbias_clean /
    mnar_sd_clean /
    mnar_rmse_clean
) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.key.size = unit(0.65, "cm"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.box.margin = margin(t = 5)
  )

mnarplot_clean