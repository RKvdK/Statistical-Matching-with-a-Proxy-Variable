# This script contains the helper functions that are required for the analysis
# This script also contains helper functions to create the plots

# ==============================================================================
# asFactorLevels()
# ==============================================================================

# Helper function to ensure the correct factor level structure
# Variable x is ensured to have the factor levels, specified in lvls

asFactorLevels <- function(x, lvls) { 
  
  # If x is numeric or integer
  # The numeric values are treated as integer indices to retrieve the corresponding factor levels from lvls
  
  if (is.numeric(x) || is.integer(x)) { 
    
    # Convert numeric values to integer indices
    # Retrieve the corresponding factor levels from lvls
    
    x <- lvls[as.integer(x)] 
    
  } 
  
  # If x is not numeric or integer
  # Treated as factor variables and the specified factor levels are imposed
  
  factor(x, levels = lvls) 
}

# ==============================================================================
# genZ()
# ==============================================================================

# Helper function to draw latent variable Z given Y and P
# To generate Z, transition matrices of Z given P and Z given Y are constructed
# These two matrices are integrated using the cia parameter to model violations of the conditional independence assumption

genZ <- function(Ypop, Ppop, Ylvl, Zlvl, Plvl, 
                   diagprob = 0.4, # Diagonal transition probability for Z given P and Z given Y
                   cia = 0) { # If zero, the conditional independence assumption holds
  
  nZ <- length(Zlvl) # Number of categories of latent variable Z
  nP <- length(Plvl) # Number of categories of proxy variable P
  
  # Construct transition matrix of Z given P
  
  # First, a transition matrix is created where all off-diagonal elements are equal
  # This equal value corresponds to (1 - diagprob) / (nP - 1) = (1-0.4) / (3-1) = 0.3
   
  tranmatZP <- matrix((1 - diagprob) / (nP - 1), 
                      nrow = nZ, ncol = nP, 
                      dimnames = list(Zlvl, Plvl)) 
  
  # Second, the diagonal elements of the transition matrix are replaced by diagprob to model the dependency of Z on P
  
  for (i in seq_len(min(nZ, nP))) { 
    
    tranmatZP[i, i] <- diagprob # Insert diagprob on the diagonal 
    
  }
  
  # Third, the transition matrix is normalized to ensure that the column sums equal 1
  
  tranmatZP <- sweep(tranmatZP, 2, colSums(tranmatZP), "/") # Normalize the transition matrix
  
  # Construct transition matrix of Z given Y
  
  nY <- length(Ylvl) # Number of categories of variable Y
  
  # First, a transition matrix is created where all off-diagonal elements are equal
  # This equal value again corresponds to (1 - diagprob) / (nY - 1) = (1-0.4) / (3-1) = 0.3
  
  tranmatZY <- matrix((1 - diagprob) / (nY - 1), nrow = nZ, ncol = nY, # Transition matrix Z given Y
                      dimnames = list(Zlvl, Ylvl))
  
  
  # Second, the diagonal elements of the transition matrix are replaced by diagprob to model the dependency of Z on Y
  
  for (i in seq_len(min(nZ, nY))) {
    
    tranmatZY[i, i] <- diagprob # Insert diagprob on the diagonal
    
  }
  
  # Third, the transition matrix is normalized to ensure that the column sums equal 1
  
  tranmatZY <- sweep(tranmatZY, 2, colSums(tranmatZY), "/") # Normalize the transition matrix
  
  # Column indices are generated to retrieve the correct columns from the transition matrices based on the values of Ppop and Ypop
  
  idxP <- if (is.numeric(Ppop) || is.integer(Ppop)) as.integer(Ppop) else as.integer(factor(Ppop, levels = Plvl))   # Column index for P
  idxY <- if (is.numeric(Ypop) || is.integer(Ypop)) as.integer(Ypop) else as.integer(factor(Ypop, levels = Ylvl)) # Column index for Y
  
  # Retrieve the corresponding conditional probability vectors for each observation in the population from the transition matrices
  
  ZPXmat <- tranmatZP[, idxP, drop = FALSE] # Retrieve the conditional probability vectors of Z for each observation in the population given P
  ZYmat  <- tranmatZY[, idxY, drop = FALSE] # Retrieve the conditional probability vectors of Z for each observation in the population given Y
  
  # To induce violations of the conditional independence assumption, the two transition matrices are integrated using the cia parameter
  
  Zmat <- (1 - cia) * ZPXmat + cia * ZYmat # Integrate cia parameter 
  Zmat <- sweep(Zmat, 2, colSums(Zmat), "/") # Normalize the resulting matrix
  
  # Based on the resulting conditional probability vectors, Z is drawn for each observation in the population using a multinomial distribution
  # p is the conditional probability vector for each observation in the population 
  # rmultinom is used to draw from a multinomial distribution with one trial 
  # which == 1 is used to retrieve the index of the category of Z that is drawn for each observation in the population
  
  Zdraw <- apply(Zmat, 2, function(p) which(rmultinom(1, 1, p) == 1)) # Draw Z given Y and P
  Zpop <- factor(Zlvl[Zdraw], levels = Zlvl) # Ensure correct factor structure
  
  # Return factor variable Z
  
  factor(Zpop, levels = Zlvl) 
  
}

# ==============================================================================
# randoff()
# ==============================================================================

# Helper function to create random off-diagonal elements for the transition matrix

randoff <- function(r) { # Function to create random off-diagonal elements
  
  w <- rexp(r) # Draw randomly values from a gamma(1,1) distribution
  w / sum(w) # Normalize such that the elements sum to 1
  
  # In this case, a gamma(1,1) distribution (default) is used to ensure positive values
  # a gamma(1,1) distribution is shape wise comparable to an uniform distribution
  # A gamma distribution is used instead of an uniform distribution to avoid zero values and is more flexible
  
}

# ==============================================================================
# marSM()
# ==============================================================================

# Helper function to compute the marginal probability of variable VAR 

# First the variable is ensured to be a factor variable 
# If factor levels are specified, these are imposed as well

marSM <- function(VAR, levels = NULL, a = 0.5) {  
  
  if (is.null(levels)) { # If there are no factor levels specified
    
    VAR <- factor(VAR) # VAR is ensured to be a factor variable          
    
  } else { # If there are factor levels specified
    
    VAR <- asFactorLevels(VAR, levels) # VAR is ensured to be a factor variable 
    
  }
  
  # Then, the marginal probability is computed using the table()-function
  # A smoothing factor is integrated to avoid cells containing zero counts
  
  tab <- table(VAR) # Compute the marginal probabilities
  (tab + a) / sum(tab + a) # Integrate the smoothing factor and normalize
  
}

# ==============================================================================
# conSM()
# ==============================================================================

# Helper function to compute the conditional probability of a variable given a second variable
# The given argument is used to condition either on the column or row variable
# If factor levels are specified, they can be imposed on the column and row variables

conSM <- function(RVAR, CVAR, given = c("col","row"),
                  rlv = NULL, 
                  clv = NULL, 
                  a = 0.5) { 
  
  given <- match.arg(given) # Condition on either the column or row variable
  
  # In case conditioning is on the column variable
  
  if (is.null(rlv)) { # If there are no factor levels specified for the row variable
    
    RVAR <- factor(RVAR) # RVAR is ensured to be a factor variable          
    
  } else { # If there are factor levels specified for the row variable
    
    RVAR <- asFactorLevels(RVAR, rlv) # Enforce factor levels to the row variable if present        
    
  }
  
  # In case conditioning is on the row variable
  
  if (is.null(clv)) { # If there are no factor levels specified for the column variable
    
    CVAR <- factor(CVAR) # CVAR is ensured to be a factor variable
    
  } else { # If there are factor levels specified for the column variable
    
    CVAR <- asFactorLevels(CVAR, clv) # Enforce factor levels to the column variable if present
    
  }
  
  # Construct joint frequency matrix with correct dimensions and factor levels
  # A smoothing factor is integrated to avoid cells containing zero counts
  
  J <- table(RVAR, CVAR)        
  J <- (J + a) / sum(J + a) # Integrate the smoothing factor and normalize
  
  # In case the distribution is conditioned on the column variable
  
  if (given == "col") { 
    
    conprob <- sweep(J, 2, colSums(J), "/") # Condition on the column variable
    
  # In case the distribution is conditioned on the row variable  
  
  } else { 
    
    conprob <- sweep(J, 1, rowSums(J), "/") # Condition on the row variable
  }
  
  # Return the conditional probability matrix
  
  return(conprob)  
}

# ==============================================================================
# orderYZ()
# ==============================================================================

# Helper function to order the YZ matrix correctly

orderYZ <- function(M, Ylvl = NULL, Zlvl = NULL) { 
  
  if (is.null(Ylvl)) { # If no Y levels are specified
    
    # Extract Y levels from row names of M or from the table structure of M
    
    Ylvl <- if (!is.null(rownames(M))) rownames(M) else sort(unique(rownames(as.table(M)))) 
    
  }
  
  if (is.null(Zlvl)) { # If no Z levels are specified
    
    # Extract Z levels from column names of M or from the table structure of M
    
    Zlvl <- if (!is.null(colnames(M))) colnames(M) else sort(unique(colnames(as.table(M))))
    
  } 
  
  # Enforce the correct order of Y and Z levels in the matrix M
  # Return the ordered matrix
  
  M[Ylvl, Zlvl, drop = FALSE] 
  
}

# ==============================================================================
# alignP()
# ==============================================================================

# Helper function to align all proxy variable levels

alignP <- function(YconP, PconX, ZconPXlist, Plvl) {  
  
  # Construct a list containing all proxy variable level vectors
  
  Psets <- c(list(colnames(YconP), rownames(PconX)), lapply(ZconPXlist, colnames)) 
  
  # Extract the intersection of the proxy variable level
  
  P <- Reduce(intersect, Psets) 
  
  # If the intersection is empty, the default proxy factor levels are inserted
  
  if (length(P) == 0L) P <- Plvl 
  
  # Enforce the intersected proxy variable levels
  # Store the resulting matrices and the used proxy variable levels in a list
  
  list( 
    YconP  = YconP[, P, drop = FALSE],
    PconX  = PconX[P, , drop = FALSE],
    ZconPX = lapply(ZconPXlist, function(M) M[, P, drop = FALSE]),
    Pused  = P
  )
}

# ==============================================================================
# ipfSafe()
# ==============================================================================

# Helper function to perform Iterative Proportional Fitting (IPF)

ipfSafe <- function(seed, # Seed matrix, the unit overlap of samples A and B in principle
                    rowtar, coltar, # Target marginal distributions of the row and column variables
                    maxit = 5000, # Maximal number of iterations
                    tol = 1e-10, # Tolerance factor to specify the break condition
                    eps = 1e-10) { # Smoothing factor to avoid cells containing zero counts
  
  # Seed matrix
  
  M <- seed # Seed matrix
  M <- M + eps # Integrate additive smoothing factor
  M <- M / sum(M) # Normalize the resulting seed matrix
  
  # Target marginal distributions
  
  r <- as.numeric(rowtar) # Target marginal distribution of the row variable 
  c <- as.numeric(coltar) # Target marginal distribution of the column variable 
  
  # If the target marginal contains a NA, then this NA is replaced by zero
  
  r[is.na(r)] <- 0 # Replace NAs in specified target row marginal by 0
  c[is.na(c)] <- 0 # Replace NAs in specified target column marginal by 0
  
  # If all elements in the target distribution are zero, the target is assumed to be uniformly distributed

  if (sum(r) == 0) { # If all elements in the target row marginal are zero
    
    r[] <- 1 / length(r) # The target row marginal is uniformly distributed
    
  }
  
  if (sum(c) == 0) { # If all elements in the target column marginal are zero
    
    c[] <- 1 / length(c) # The target column marginal is uniformly distributed
    
  }
  
  # Initialize iteration counter 
  
  i <- 0 
  
  for (it in 1:maxit) {
    
    i <- it # Iteration count
    
    # Row rescaling
    
    rs <- rowSums(M) # Compute the (old) row totals
    rs[rs == 0] <- NA_real_ # Replace zeroes to avoid division by zero
    M <- sweep(M, 1, r / rs, "*") # Rescale rows using the row target
    
    # Column rescaling
    
    cs <- colSums(M) # Compute the (old) column totals
    cs[cs == 0] <- NA_real_  # Replace zeroes to avoid division by zero
    M <- sweep(M, 2, c / cs, "*") # Rescale columns using the column target
    
    # Compute the difference between the rescaled and target marginals
    
    row_err <- abs(rowSums(M) - r) # Compute the difference between the row totals and the target
    col_err <- abs(colSums(M) - c) # Compute the difference between the column totals and the target
    
    # Apply breaking condition based on the pre-defined tolerance factor
    
    if (max(row_err, na.rm = TRUE) < tol && max(col_err, na.rm = TRUE) < tol) { 
      break 
    }
  }
  
  # Normalize resulting joint distribution
  
  M <- M / sum(M) 
  
  # Store the resulting joint distribution and the number of iterations in a list
  
  out <- list(
    matrix = M,                
    iterations = i              
  )
  
  return(out)
}

# ==============================================================================
# iterdat()
# ==============================================================================

# Helper function to convert the RUN data such that each row represents an iteration

iterdat <- function(runres, estimator = "DRE", include_bias = TRUE) {
  
  # For every estimator, the estimated joint distribution of Y and Z is extracted
  # The true joint distribution of Y and Z is extracted by correcting the estimated joint distribution for the cell bias
  
  # Doubly robust estimator 
  
  if (estimator == "DRE") { # 
    
    theta <- runres$DRE[[1]]$theta_DRE # Extract the estimated joint distribution of Y and Z
    YZpop <- runres$DRE[[1]]$mean_theta_DRE - runres$DRE[[1]]$bias_cell_DRE # Extract the true joint distribution of Y and Z
    
  }
  
  # Iterative proportional fitting estimator
  
  if (estimator == "IPF") { 
    
    theta <- runres$IPF[[1]]$theta_IPF # Extract the estimated joint distribution of Y and Z
    YZpop <- runres$IPF[[1]]$mean_theta_IPF - runres$IPF[[1]]$bias_cell_IPF # Extract the true joint distribution of Y and Z
  }
  
  # CIA proxy estimator
  
  if (estimator == "PROXY") { 
    
    theta <- runres$PROXY[[1]]$theta_PROXY # Extract the estimated joint distribution of Y and Z
    YZpop <- runres$PROXY[[1]]$mean_theta_PROXY - runres$PROXY[[1]]$bias_cell_PROXY # Extract the true joint distribution of Y and Z
    
  }
  
  # Store the resulting data in a data frame 
  
  df <- as.data.frame.table(theta) # Convert the estimated joint distribution of Y and Z to a data frame
  names(df) <- c("iter", "Y", "Z", "theta") # Rename the columns of the data frame
  df$iter <- as.integer(df$iter) # Ensure that the iteration variable is an integer
  
  # If include_bias = TRUE, additional data is added to the data frame to compute the bias for each cell and iteration
  
  if (include_bias) {
    
    df$true <- mapply(function(y, z) YZpop[y, z], df$Y, df$Z) # Retrieve the true joint distribution of Y and Z for each cell
    df$bias <- df$theta - df$true # Compute the bias for each cell and iteration
    
  }
  
  return(df) 
  
}

# ==============================================================================
# concell()
# ==============================================================================

# Helper function to evaluate convergence on the cell level

concell <- function(condat, scenario) {
  
  conv <- condat %>% 
    
    # Group the data by cell and arrange by iteration within each cell
    
    group_by(Y, Z) %>% # Group by cell 
    arrange(iter, .by_group = TRUE) %>% # Arrange by iteration within each cell
    
    # Compute the running bias and running standard deviation for each cell and iteration
    
    mutate( 
      
      # Compute running bias for each cell
      
      running_bias = cumsum(bias) / seq_along(bias), 
     
      # Compute running standard deviation for each cell
      
      running_sd = sapply(seq_along(theta), function(i){ 
        
        # Standard deviation calculation is only performed if there are at least 2 iterations, otherwise NA is returned
        
        if (i == 1) return(NA_real_) 
        
        # Compute the running standard deviation for each cell
        
        sd(theta[1:i], na.rm = TRUE) 
      }),
      
      # Create identifier for each cell
      
      cell = paste0("(", Y, ", ", Z, ")") 
      
    ) %>%
    
    ungroup() # Ungroup the data
  
  # Create running bias plot
  
  p_bias <- ggplot(conv, aes(x = iter, y = running_bias, color = cell)) +
    geom_line() + 
    labs(
      x = "MC iterations",
      y = "Running bias",
      color = "Cell",
      title = paste("Running bias per cell -", scenario)
    ) +
    theme_bw()
  
  # Create running standard deviation plot
  
  p_sd <- ggplot(conv, aes(x = iter, y = running_sd, color = cell)) +
    geom_line() + 
    labs(
      x = "MC iterations",
      y = "Running SD",
      color = "Cell",
      title = paste("Running SD per cell -", scenario)
    ) +
    theme_bw()
  
  # Store the resulting plots and data in a list
  
  list( 
    bias_plot = p_bias,
    sd_plot = p_sd,
    conv_data = conv
  )
  
}

# ==============================================================================
# barplotSM()
# ==============================================================================

# Helper function to create bar plots comparing the performance of the estimators across scenarios

barplotSM <- function(df, cols, value_name, title) {
  
  plot_data <- df %>% 
    
    # Select the scenario and estimator columns
    
    select(scenario, all_of(cols)) %>% 
    
    # Convert into long format
    
    pivot_longer( 
      cols = -scenario,
      names_to = "estimator",
      values_to = "value"
    ) %>%
    
    mutate(
      
      estimator = gsub(".*_", "", estimator) , # Extract estimator name
      estimator = factor(estimator, levels = c("DRE", "IPF", "PROXY")), # Convert estimator into a factor variable
      scenario = factor(scenario, levels = unique(df$scenario)) # Convert scenario into a factor variable
    )
  
  # Identify the best performing estimator in each scenario
  
  winners <- plot_data %>% 
    
    # Group by scenario
    
    group_by(scenario) %>% 
    
    # Identify the row with the minimum value for each scenario
    
    slice_min(order_by = value, n = 1, with_ties = FALSE) %>% 
    
    # Ungroup the data
    
    ungroup() %>% 
    
    # Adjust position for the asterisk based on the estimator
    
    mutate(
      x_pos = case_when( 
        estimator == "DRE" ~ as.numeric(scenario) - 0.3,
        estimator == "IPF" ~ as.numeric(scenario),
        estimator == "PROXY" ~ as.numeric(scenario) + 0.3
      )
    )
  
  # Create the bar plot
  
  ggplot(plot_data, aes(x = scenario, y = value, fill = estimator, group = estimator)) +
    
    # Bar plot visualizing the performance of the estimators across scenarios
    
    geom_col(
      position = position_dodge(width = 0.9),
      width = 0.8,
      color = "black",
      linewidth = 0.35
    ) +
    
    # Add asterisks to indicate the best performing estimator in each scenario
    
    geom_text( 
      data = winners,
      aes(x = x_pos, y = value, label = "*"),
      vjust = -0.4, # Adjust vertical position of the asterisk
      size = 6, # Size of the asterisk
      inherit.aes = TRUE # Inherit aesthetics from the main plot
      ) +
    labs(
      x = "Scenario",
      y = value_name,
      fill = "Estimator",
      title = title
    ) +
    theme_bw()
}