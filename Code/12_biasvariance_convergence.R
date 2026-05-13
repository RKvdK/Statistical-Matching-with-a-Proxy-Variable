# This script evaluates the bias and variance convergence
# As function of the number of Monte Carlo iterations

# ==============================================================================
# Benchmark Bias Estimations
# ==============================================================================

# Compute Benchmark Values

# A1_MC2000 <- RUN(A1, 2000, keep_iter = TRUE)
# D7_MC2000 <- RUN(D7, 2000, keep_iter = TRUE)
# H16_MC2000 <- RUN(H16, 2000, keep_iter = TRUE)

# Compute iteration-level convergence data

# A1_condat <- iterdat(A1_MC2000)
# D7_condat <- iterdat(D7_MC2000)
# H16_condat <- iterdat(H16_MC2000)

# Save iteration-convergence data

# saveRDS(A1_condat, "A1_benchmark")
# saveRDS(D7_condat, "D7_benchmark")
# saveRDS(H16_condat, "H16_benchmark")

# Load iteration-convergence data

A1_condat <- readRDS("A1_benchmark")
D7_condat <- readRDS("D7_benchmark")
H16_condat <- readRDS("H16_benchmark")

# ==============================================================================
# Convergence on the cell level
# ==============================================================================

# Compute running bias and running standard deviation for every Y-Z cell

A1_concell <- concell(A1_condat, "A1")
D7_concell <- concell(D7_condat, "D7")
H16_concell <- concell(H16_condat, "H16")

# Visualize convergence of running bias

A1_concell$bias_plot
D7_concell$bias_plot
H16_concell$bias_plot

# Visualize convergence of running standard deviation

A1_concell$sd_plot
D7_concell$sd_plot
H16_concell$sd_plot

# Modify bias plots

p_A1_bias <- A1_concell$bias_plot +
  scale_color_viridis_d( # Colorblind-friendly palette
    option = "C", 
    begin = 0.1,
    end = 0.9
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank()
  )

p_D7_bias <- D7_concell$bias_plot +
  scale_color_viridis_d( # Colorblind-friendly palette
    option = "C",
    begin = 0.1,
    end = 0.9
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank()
  )

p_H16_bias <- H16_concell$bias_plot +
  scale_color_viridis_d( # Colorblind-friendly palette
    option = "C",
    begin = 0.1,
    end = 0.9
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

# Modify sd plots

p_A1_sd <- A1_concell$sd_plot +
  scale_color_viridis_d( # Colorblind-friendly palette
    option = "C",
    begin = 0.1,
    end = 0.9
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  scale_y_continuous(position = "right") +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank()
  )

p_D7_sd <- D7_concell$sd_plot +
  scale_color_viridis_d( # Colorblind-friendly palette
    option = "C",
    begin = 0.1,
    end = 0.9
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  scale_y_continuous(position = "right") +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank()
  )

p_H16_sd <- H16_concell$sd_plot +
  scale_color_viridis_d( # Colorblind-friendly palette
    option = "C",
    begin = 0.1,
    end = 0.9
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  scale_y_continuous(position = "right") +
  theme(
    legend.position = "bottom",
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

# Column labels

col_bias <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Bias", size = 5, fontface = "plain") +
  theme_void()

col_sd <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Standard deviation", size = 5, fontface = "plain") +
  theme_void()

# Row labels

row_A1 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "A1", size = 5, fontface = "plain") +
  theme_void()

row_D7 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "D7", size = 5, fontface = "plain") +
  theme_void()

row_H16 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "H16", size = 5, fontface = "plain") +
  theme_void()

# Header and rows

left_width <- 0.14
main_width <- 1

header_row <- plot_spacer() + col_bias + col_sd +
  plot_layout(widths = c(left_width, main_width, main_width))

r1 <- row_A1 + p_A1_bias + p_A1_sd +
  plot_layout(widths = c(left_width, main_width, main_width))

r2 <- row_D7 + p_D7_bias + p_D7_sd +
  plot_layout(widths = c(left_width, main_width, main_width))

r3 <- row_H16 + p_H16_bias + p_H16_sd +
  plot_layout(widths = c(left_width, main_width, main_width))


x_axis_label <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Number of Monte Carlo draws", size = 6) +
  theme_void()

# Combined figure

final_plot <- (header_row / r1 / r2 / r3) +
  plot_layout(
    heights = c(0.08, 1, 1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

final_plot
