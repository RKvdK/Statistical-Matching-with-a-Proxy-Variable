# In this script, all specified functions and objects are combined to retrieve the actual results
# The SIM, MC and RUN functions are required for this script to run adequately

# ==============================================================================
# Simulate and save the data
# ==============================================================================

# Run the simulation with X Monte Carlo iterations per parameter combination

# MC50 <- RUN(simgrid, MCnum = 50) 

# Save the simulation results

# saveRDS(MC50, file = "MC50.rds") 

# ==============================================================================
# Load the simulated data
# ==============================================================================

# Load the simulation results

MC50 <- readRDS("MC50.rds")