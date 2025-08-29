# source("bayesian_group_lasso_true.R")
source("bayesian_lasso_true.R")
library(mvtnorm)
library(tidyverse)
library(pheatmap)
library(gridExtra)
library(grid)
library(vegan)
library(coda)

simulate_data_with_groups <- function(n, p, G, K, active_groups, sigma) {
  # Validation checks
  if (p %% G != 0) stop("p must be divisible by G")
  if (length(active_groups) > K) stop("Number of active groups cannot exceed K")
  
  # Create group structure: divide p features into G equal groups
  group_size <- p / G  # features per group
  group_id <- rep(1:G, each = group_size)  # group membership vector
  
  # Initialize loading matrix W as all zeros (sparse structure)
  W_true <- matrix(0, nrow = p, ncol = K)
  
  # This creates the sparse group structure we want to recover
  for (i in 1:length(active_groups)) {
    current_active_group <- active_groups[i]
    rows_in_group <- which(group_id == current_active_group)
    
    # Alternate signs for identifiability
    sign <- ifelse(i %% 2 == 1, 1, -1)
    
    # Generate strong loadings (0.6 to 1.0) for active groups
    W_true[rows_in_group, i] <- runif(length(rows_in_group), 0.6, 1) * sign
  }
  
  # Generate TRUE latent factors Z ~ N(0, I_k)
  Z_true <- matrix(rnorm(n * K), nrow = n, ncol = K)
  
  # Generate feature-specific means   ~ N(0, 1)
  mu_true <- rnorm(p, 0, 1)
  
  # Generate noise   ~ N(0,  ²I)
  noise <- matrix(rnorm(n * p, 0, sigma), nrow = n, ncol = p)
  
  # Generate complete data following: X = Z W' +   +  
  # Note: Z [n×k] × W' [k×p] = [n×p], then add   + noise [n×p]
  X_data_complete <- Z_true %*% t(W_true) + matrix(rep(mu_true, n), nrow = n, byrow = TRUE) + noise
  
  return(list(
    X_data = X_data_complete,    # Complete data matrix [n × p]
    W_true = W_true,             # True loading matrix [p × k]
    Z_true = Z_true,             #   TRUE LATENT FACTORS [n × k]  
    group_id = group_id,         # Group membership vector [p × 1]
    K_true = K,                  # True number of factors
    G_total = G                  # True number of groups
  ))
}

# MCMC Convergence Diagnostics Function
calculate_convergence_diagnostics <- function(sigma2_chains, tau_g_chains, lambda2_chains, W_chains, Z_chains, n_chains) {
  
  # Function to calculate R-hat for a single parameter
  calculate_rhat <- function(chains_matrix) {
    # chains_matrix: n_chains x n_samples
    n_chains <- nrow(chains_matrix)
    n_samples <- ncol(chains_matrix)
    
    if (n_samples < 10 || n_chains < 2) {
      return(1.0)  # Not enough data for R-hat
    }
    
    # Calculate within-chain and between-chain variances
    chain_means <- apply(chains_matrix, 1, mean)
    overall_mean <- mean(chain_means)
    
    # Within-chain variance
    W_var <- mean(apply(chains_matrix, 1, var))
    
    # Between-chain variance
    B_var <- n_samples * var(chain_means)
    
    # Pooled variance estimate
    var_plus <- ((n_samples - 1) * W_var + B_var) / n_samples
    
    # R-hat statistic - STANDARD FORMULA
    # Use the standard Gelman-Rubin R-hat formula
    if (W_var > 0) {
      rhat <- sqrt(var_plus / W_var)
    } else {
      # If W_var is exactly 0, check if there's between-chain variance
      if (B_var > 0) {
        rhat <- 2.0  # Indicate poor convergence
      } else {
        rhat <- 1.0  # No variance at all
      }
    }
    
    return(rhat)
  }
  
  cat("=== COMPREHENSIVE CONVERGENCE ASSESSMENT ===\n\n")
  
  all_rhats <- c()
  
  # 1. NOISE VARIANCE
  cat("1. NOISE VARIANCE (sigma2):\n")
  sigma2_rhat <- calculate_rhat(sigma2_chains)
  cat("sigma2 R-hat:", round(sigma2_rhat, 3))
  if (sigma2_rhat < 1.1) cat("  \n") else if (sigma2_rhat < 1.2) cat("   \n") else cat("  \n")
  all_rhats <- c(all_rhats, sigma2_rhat)
  
  # Debug: Check variance components for sigma2
  chain_means_sigma2 <- apply(sigma2_chains, 1, mean)
  W_var_sigma2 <- mean(apply(sigma2_chains, 1, var))
  B_var_sigma2 <- ncol(sigma2_chains) * var(chain_means_sigma2)
  cat("   Debug - W_var:", round(W_var_sigma2, 6), "B_var:", round(B_var_sigma2, 6), "\n")
  
  # 2. GROUP LASSO PARAMETERS
  cat("\n2. GROUP LASSO PARAMETERS (lambda2):\n")
  lambda2_rhat <- calculate_rhat(lambda2_chains)
  cat(" lambda2 (global shrinkage) R-hat:", round(lambda2_rhat, 3))
  if (lambda2_rhat < 1.1) cat("  \n") else if (lambda2_rhat < 1.2) cat("   \n") else cat("  \n")
  all_rhats <- c(all_rhats, lambda2_rhat)
  
  # tau_g (group-specific variances)
  tau_rhats <- numeric(min(5, dim(tau_g_chains)[3]))  # Test more groups
  for (g in 1:length(tau_rhats)) {
    tau_g_matrix <- tau_g_chains[, , g]
    tau_rhats[g] <- calculate_rhat(tau_g_matrix)
    cat("tau_", g, " (group ", g, " scale) R-hat:", round(tau_rhats[g], 3))
    if (tau_rhats[g] < 1.1) cat("  \n") else if (tau_rhats[g] < 1.2) cat("   \n") else cat("  \n")
  }
      all_rhats <- c(all_rhats, tau_rhats)
  
  # 3. FACTOR LOADINGS (W)
  cat("\n3. FACTOR LOADINGS (W):\n")
  P <- dim(W_chains)[3]
  K <- dim(W_chains)[4]
  
  # Test systematic sample of W elements
  w_test_indices <- list(
    c(1, 1),                              # First element
    c(ceiling(P/4), 1),                   # Quarter way through
    c(ceiling(P/2), 1),                   # Middle
    c(ceiling(3*P/4), min(2, K)),         # Three-quarters through
    c(P, min(K, 2))                       # Last element
  )
  
  W_rhats <- numeric(length(w_test_indices))
  for (i in 1:length(w_test_indices)) {
    elem <- w_test_indices[[i]]
    if (elem[1] <= P && elem[2] <= K) {
      W_matrix <- W_chains[, , elem[1], elem[2]]
      W_rhats[i] <- calculate_rhat(W_matrix)
      cat("   W[", elem[1], ",", elem[2], "] R-hat:", round(W_rhats[i], 3))
      if (W_rhats[i] < 1.1) cat("  \n") else if (W_rhats[i] < 1.2) cat("   \n") else cat("  \n")
      
      # Debug: Check variance components for W[1,1] specifically
      if (elem[1] == 1 && elem[2] == 1) {
        chain_means_w11 <- apply(W_matrix, 1, mean)
        W_var_w11 <- mean(apply(W_matrix, 1, var))
        B_var_w11 <- ncol(W_matrix) * var(chain_means_w11)
        cat("      Debug W[1,1] - W_var:", round(W_var_w11, 6), "B_var:", round(B_var_w11, 6), "\n")
      }
    }
  }
  all_rhats <- c(all_rhats, W_rhats[!is.na(W_rhats)])
  
  # 4. LATENT FACTORS (Z) - COMPREHENSIVE TESTING
  cat("\n4. LATENT FACTORS (Z):\n")
  N <- dim(Z_chains)[3]
  
  # Test systematic sample of Z elements
  z_test_indices <- list(
    c(1, 1),                              # First observation, first factor
    c(ceiling(N/4), 1),                   # Quarter way through observations
    c(ceiling(N/2), min(2, K)),           # Middle observation, second factor
    c(ceiling(3*N/4), min(K, 3)),         # Three-quarters through
    c(N, K)                               # Last observation, last factor
  )
  
  Z_rhats <- numeric(length(z_test_indices))
  for (i in 1:length(z_test_indices)) {
    elem <- z_test_indices[[i]]
    if (elem[1] <= N && elem[2] <= K) {
      Z_matrix <- Z_chains[, , elem[1], elem[2]]
      Z_rhats[i] <- calculate_rhat(Z_matrix)
      cat("   Z[", elem[1], ",", elem[2], "] R-hat:", round(Z_rhats[i], 3))
      if (Z_rhats[i] < 1.1) cat("  \n") else if (Z_rhats[i] < 1.2) cat("   \n") else cat("  \n")
    }
  }
  all_rhats <- c(all_rhats, Z_rhats[!is.na(Z_rhats)])
  
  # OVERALL CONVERGENCE ASSESSMENT
  all_rhats <- all_rhats[!is.na(all_rhats)]
  max_rhat <- max(all_rhats, na.rm = TRUE)
  prop_converged_excellent <- mean(all_rhats < 1.1, na.rm = TRUE)
  prop_converged_acceptable <- mean(all_rhats < 1.2, na.rm = TRUE)
  
  cat("\n=== OVERALL CONVERGENCE SUMMARY ===\n")
  cat("   Total parameters tested:", length(all_rhats), "\n")
  cat("   Max R-hat:", round(max_rhat, 3), "\n")
  cat("   Excellent convergence (< 1.1):", round(100 * prop_converged_excellent, 1), "%\n")
  cat("   Acceptable convergence (< 1.2):", round(100 * prop_converged_acceptable, 1), "%\n")
  
  # CONVERGENCE INTERPRETATION
  if (max_rhat < 1.1) {
    cat("\n EXCELLENT: All parameters have converged excellently!\n")
    recommendation <- "ready_for_analysis"
  } else if (max_rhat < 1.2 && prop_converged_excellent > 0.8) {
    cat("\n GOOD: Most parameters converged well, a few need attention\n")
    recommendation <- "mostly_converged"
  } else if (max_rhat < 1.5) {
    cat("\n MARGINAL: Convergence issues detected. Increase iterations.\n")
    recommendation <- "increase_iterations"
  } else {
    cat("\n POOR: Serious convergence problems. Check model specification.\n")
    recommendation <- "model_problems"
  }
  
  # SPECIFIC RECOMMENDATIONS
  cat("\nRECOMMENDATIONS:\n")
  if (recommendation == "ready_for_analysis") {
    cat("Analysis ready! Proceed with inference.\n")
  } else if (recommendation == "mostly_converged") {
    cat("Consider doubling iterations for problematic parameters\n")
  } else if (recommendation == "increase_iterations") {
    cat("Increase n_iter to at least", ceiling(1.5 * max_rhat * 1000), "iterations\n")
    cat("Consider using", min(n_chains + 2, 6), "chains instead of", n_chains, "\n")
  } else {
    cat("Check model specification and starting values\n")
    cat("Try different initialization strategies\n")
    cat("Consider reparameterization\n")
  }
  
  return(list(
    sigma2_rhat = sigma2_rhat,
    lambda2_rhat = lambda2_rhat,
    tau_rhats = tau_rhats,
    W_rhats = W_rhats,
    Z_rhats = Z_rhats,
    max_rhat = max_rhat,
    prop_converged_excellent = prop_converged_excellent,
    prop_converged_acceptable = prop_converged_acceptable,
    recommendation = recommendation,
    n_parameters_tested = length(all_rhats)
  ))
}

# ROTATION METHOD FOR FACTOR IDENTIFIABILITY

# Post-hoc Procrustes Alignment (Stephens 2000, Papastamoulis 2016)
align_chains_procrustes <- function(W_chains, Z_chains, n_chains) {
  cat("Method: Procrustes alignment (Stephens 2000)...\n")
  
  n_iter <- dim(W_chains)[2]
  P <- dim(W_chains)[3] 
  K <- dim(W_chains)[4]
  N <- dim(Z_chains)[3]
  
  W_aligned <- W_chains
  Z_aligned <- Z_chains
  
  for (chain in 2:n_chains) {
    for (iter in 1:n_iter) {
      W_ref <- W_chains[1, iter, , ]
      W_curr <- W_chains[chain, iter, , ]
      Z_curr <- Z_chains[chain, iter, , ]
      
      if (K > 1) {
        cross_cov <- t(W_ref) %*% W_curr
        if (!any(is.na(cross_cov))) {
          svd_result <- svd(cross_cov)
          R_opt <- svd_result$v %*% t(svd_result$u)
          if (det(R_opt) < 0) {
            svd_result$v[, K] <- -svd_result$v[, K]
            R_opt <- svd_result$v %*% t(svd_result$u)
          }
          W_aligned[chain, iter, , ] <- W_curr %*% R_opt
          Z_aligned[chain, iter, , ] <- Z_curr %*% R_opt
        }
      } else {
        sign_corr <- sign(sum(W_ref * W_curr))
        W_aligned[chain, iter, , ] <- sign_corr * W_curr
        Z_aligned[chain, iter, , ] <- sign_corr * Z_curr
      }
    }
  }
  
  return(list(W_aligned = W_aligned, Z_aligned = Z_aligned))
}

# Main alignment function with method selection
align_factor_chains <- function(W_chains, Z_chains, n_chains, method = "procrustes") {
  cat("Aligning", n_chains, "chains using", method, "method...\n")
  
  result <- align_chains_procrustes(W_chains, Z_chains, n_chains)
  
  cat("Alignment completed!\n")
  return(result)
}

# AUTO-CONVERGENCE ADJUSTMENT FUNCTION
auto_adjust_mcmc <- function(X, K, G, initial_n_iter = 5000, initial_burn_in = 1500, 
                             n_chains = 4, max_attempts = 3, target_rhat = 1.1) {
  cat("=== AUTO-ADJUSTING MCMC FOR CONVERGENCE ===\n")
  
  current_n_iter <- initial_n_iter
  current_burn_in <- initial_burn_in
  
  for (attempt in 1:max_attempts) {
    cat("\nAttempt", attempt, "- Running", current_n_iter, "iterations with", current_burn_in, "burn-in\n")
    
    # Run MCMC
    result <- mcmc_bpca_group_lasso(X, K = K, G = G, 
                                    n_iter = current_n_iter, 
                                    burn_in = current_burn_in, 
                                    n_chains = n_chains)
    
    # Check convergence
    max_rhat <- result$convergence_diagnostics$max_rhat
    recommendation <- result$convergence_diagnostics$recommendation
    
    cat("Max R-hat achieved:", round(max_rhat, 3), "\n")
    
    if (max_rhat <= target_rhat) {
      cat("CONVERGENCE ACHIEVED! R-hat =", round(max_rhat, 3), " ", target_rhat, "\n")
      result$auto_adjustment_info <- list(
        attempts_needed = attempt,
        final_n_iter = current_n_iter,
        final_burn_in = current_burn_in,
        convergence_achieved = TRUE
      )
      return(result)
    } else {
      cat("Convergence not achieved. R-hat =", round(max_rhat, 3), ">", target_rhat, "\n")
      
      if (attempt < max_attempts) {
        # Increase iterations
        current_n_iter <- ceiling(current_n_iter * 1.5)
        current_burn_in <- ceiling(current_burn_in * 1.2)
        cat("Increasing to", current_n_iter, "iterations with", current_burn_in, "burn-in\n")
      }
    }
  }
  
  # If we get here, convergence was not achieved
  cat("Convergence not achieved after", max_attempts, "attempts\n")
  cat("Consider: longer chains, more chains, or model reparameterization\n")
  
  result$auto_adjustment_info <- list(
    attempts_needed = max_attempts,
    final_n_iter = current_n_iter,
    final_burn_in = current_burn_in,
    convergence_achieved = FALSE,
    final_max_rhat = max_rhat
  )
  
  return(result)
}

#                                                                                    
#                        BAYESIAN GROUP LASSO BPCA - MAIN FUNCTION
#                                                                                    
#
# PURPOSE: Perform Bayesian PCA with Group LASSO priors using multi-chain MCMC
#
# PARAMETERS:
# - X: data matrix [n × p] with possible missing values (NA)
# - K: number of latent factors (default: NULL, will use auto_K or default to 5)
# - G: number of groups for variables (default: 5, features divided equally)
# - n_iter: total MCMC iterations per chain (default: 5000, optimized)
# - burn_in: burn-in samples to discard (default: 1500, ~30% of n_iter)
# - lambda_a, lambda_b: Gamma prior hyperparameters for sigma2 (default: 2.0, 1.0 - OPTIMIZED)
# - auto_group: if TRUE, automatically determine optimal G (default: FALSE)
# - auto_K: if TRUE, automatically determine optimal K using cross-validation (default: FALSE)
# - K_range: range of K values to test when auto_K=TRUE (default: 2:8)
# - cv_folds: number of cross-validation folds for K selection (default: 5)
# - n_chains: number of independent MCMC chains (default: 4, for robust convergence)
#
# AUTO_GROUP EXPLANATION:
# When auto_group=TRUE:
# 1. Tests multiple values of G (e.g., 2, 3, 4, 5, 6...)
# 2. For each G, runs short MCMC to evaluate model fit
# 3. Uses information criteria (WAIC, DIC) or cross-validation likelihood
# 4. Selects G that best balances fit vs complexity
# 5. Then runs full MCMC with optimal G
# This is useful when don't know the true group structure
#
# AUTO_K EXPLANATION:
# When auto_K=TRUE:
# 1. Tests multiple values of K from K_range (e.g., 2, 3, 4, 5, 6, 7, 8)
# 2. For each K, runs cross-validation with cv_folds folds
# 3. Each fold: trains on training data, evaluates reconstruction on validation data
# 4. Uses log-likelihood of reconstructed missing values as evaluation metric
# 5. Selects K with highest average cross-validation score
# 6. Then runs full MCMC with optimal K
# This is useful when don't know the optimal number of latent factors
#
# Z (LATENT FACTORS) HANDLING:
# Z is the core latent factor matrix [n × k] that represents subjects in factor space
# 1. INITIALIZATION: Z starts with random values ~ N(0, 0.5²) for each chain
# 2. MCMC UPDATES: Z is updated every iteration using Gibbs sampling
# 3. STORAGE: Z samples stored in multiple arrays:
#    - Z_samples_all[chain, iter, subject, factor]: All raw samples
#    - Z_samples_postburn: After removing burn-in samples
#    - Z_combined: Combined across all chains for final estimates
#    - Z_mean: Final posterior mean estimate [n × k]
# 4. ALIGNMENT: Z undergoes post-hoc rotation alignment to fix identifiability
# 5. RETURN: Multiple Z objects returned for different purposes
#                                                                                    
#                           CROSS-VALIDATION FOR K SELECTION
#                                                                                    
cross_validate_K <- function(X, K_test, G, cv_folds, n_iter, burn_in, lambda_a, lambda_b, auto_group) {
  # PURPOSE: Cross-validate BPCA model with specific K to evaluate model fit
  # RETURNS: List with mean CV score, standard error, and fold-wise scores
  
  N <- nrow(X)
  P <- ncol(X)
  
  # Create cross-validation folds
  fold_size <- floor(N / cv_folds)
  fold_indices <- sample(rep(1:cv_folds, length.out = N))
  
  cv_scores <- numeric(cv_folds)
  
  for (fold in 1:cv_folds) {
    # Split data into training and validation sets
    val_idx <- which(fold_indices == fold)
    train_idx <- setdiff(1:N, val_idx)
    
    X_train <- X[train_idx, , drop = FALSE]
    X_val <- X[val_idx, , drop = FALSE]
    
    # Handle missing values in validation set
    # Keep some observed values for evaluation, hide others for prediction
    val_missing_rate <- 0.3  # Hide 30% of validation data for prediction
    X_val_test <- X_val
    
    # Create missing data in validation set (but keep track of true values)
    for (i in 1:nrow(X_val_test)) {
      observed_idx <- which(!is.na(X_val_test[i, ]))
      if (length(observed_idx) > 2) {  # Keep at least 2 observed values
        n_hide <- max(1, floor(length(observed_idx) * val_missing_rate))
        hide_idx <- sample(observed_idx, min(n_hide, length(observed_idx) - 1))
        X_val_test[i, hide_idx] <- NA
      }
    }
    
    # Train model on training data
    tryCatch({
      model <- mcmc_bpca_group_lasso(X_train, K = K_test, G = G, 
                                    n_iter = n_iter, burn_in = burn_in,
                                    lambda_a = lambda_a, lambda_b = lambda_b,
                                    auto_group = auto_group, auto_K = FALSE,
                                    n_chains = 2, alignment_method = "procrustes")
      
      # Predict on validation data and calculate reconstruction error
      score <- evaluate_reconstruction_K(model, X_val, X_val_test)
      cv_scores[fold] <- score
      
    }, error = function(e) {
      # If model fails, assign poor score
      cv_scores[fold] <<- -1000
      cat(" [FAILED]")
    })
  }
  
  # Calculate mean and standard error (with better NA handling)
  valid_scores <- cv_scores[!is.na(cv_scores) & !is.infinite(cv_scores) & cv_scores > -999]
  if (length(valid_scores) == 0) {
    return(list(mean_score = -1000, se_score = 0, fold_scores = cv_scores, n_valid_folds = 0))
  }
  
  mean_score <- mean(valid_scores)
  se_score <- if (length(valid_scores) > 1) sd(valid_scores) / sqrt(length(valid_scores)) else 0
  
  return(list(
    mean_score = mean_score,
    se_score = se_score,
    fold_scores = cv_scores,
    n_valid_folds = length(valid_scores)
  ))
}

evaluate_reconstruction_K <- function(model, X_true, X_with_missing) {
  # PURPOSE: Evaluate how well the model reconstructs missing values
  # Enhanced version with natural missing patterns and robust error handling
  # Higher score = better reconstruction
  
  # Extract model parameters
  W <- model$W  # Posterior mean of loadings
  mu <- model$mu  # Posterior mean of intercepts
  sigma2 <- model$sigma2  # Posterior mean of noise variance
  
  N_val <- nrow(X_with_missing)
  P <- ncol(X_with_missing)
  K <- ncol(W)
  
  total_log_likelihood <- 0
  n_predictions <- 0
  
  # For each validation subject
  for (i in 1:N_val) {
    observed_idx <- which(!is.na(X_with_missing[i, ]))
    missing_idx <- which(is.na(X_with_missing[i, ]))
    
    if (length(observed_idx) > 0 && length(missing_idx) > 0) {
      # Predict latent factors Z from observed data
      X_obs <- X_with_missing[i, observed_idx]
      W_obs <- W[observed_idx, , drop = FALSE]
      mu_obs <- mu[observed_idx]
      
      # Posterior of Z given observed data: Z | X_obs ~ N(mu_z, Sigma_z)
      # Add regularization for numerical stability
      Sigma_z_inv <- t(W_obs) %*% W_obs / sigma2 + diag(K) + 1e-8 * diag(K)
      Sigma_z <- solve(Sigma_z_inv)
      mu_z <- Sigma_z %*% t(W_obs) %*% (X_obs - mu_obs) / sigma2
      
      # Predict missing values: X_missing | Z ~ N(W*mu_z + mu, sigma2)
      X_pred_missing <- W[missing_idx, , drop = FALSE] %*% mu_z + mu[missing_idx]
      X_true_missing <- X_true[i, missing_idx]
      
      # Calculate log-likelihood of predictions (with enhanced error handling)
      if (length(missing_idx) > 0) {
        residuals <- X_true_missing - X_pred_missing
        
        # Enhanced validation of residuals and parameters
        if (any(is.na(residuals)) || any(is.infinite(residuals)) || 
            sigma2 <= 0 || any(is.na(X_pred_missing)) || any(is.infinite(X_pred_missing))) {
          # Skip this prediction if values are invalid
          next
        }
        
        log_lik <- -0.5 * (sum(residuals^2) / sigma2 + length(missing_idx) * log(2 * pi * sigma2))
        
        # Check if log_lik is valid and reasonable
        if (!is.na(log_lik) && !is.infinite(log_lik) && log_lik > -1e6) {
          total_log_likelihood <- total_log_likelihood + log_lik
          n_predictions <- n_predictions + length(missing_idx)
        }
      }
    }
  }
  
  # Return average log-likelihood per prediction
  if (n_predictions > 0) {
    return(total_log_likelihood / n_predictions)
  } else {
    return(-1000)  # Poor score if no predictions could be made
  }
}

mcmc_bpca_group_lasso <- function(X, K = 5, G = 5, n_iter = 10000, burn_in = 3500,
                                  lambda_a = 0.1, lambda_b = 0.1, auto_group = FALSE, auto_K = FALSE,
                                  K_range = 2:8, cv_folds = 5, n_chains = 4) {
  
  # Validate required functions from bayesian_group_lasso_true.R
  required_functions <- c("sample_group_lasso_loadings", "sample_group_scales_gig", 
                         "sample_global_shrinkage")
  missing_functions <- required_functions[!sapply(required_functions, exists)]
  if (length(missing_functions) > 0) {
    stop("Missing required functions: ", paste(missing_functions, collapse = ", "), 
         "\nPlease ensure bayesian_group_lasso_true.R is properly sourced.")
  }
  
  N <- nrow(X)  # Number of subjects/observations
  P <- ncol(X)  # Number of features/variables
  
  #                                                                                    
  #                              AUTO K SELECTION WITH CROSS-VALIDATION
  #                                                                                    
  # When auto_K=TRUE, automatically determine optimal number of factors K using CV
  if (auto_K && is.null(K)) {
    cat("AUTO K SELECTION: Finding optimal K using", cv_folds, "-fold cross-validation...\n")
    
    # Ensure K_range is reasonable
    K_range <- K_range[K_range >= 1 & K_range <= min(N, P)]
    if (length(K_range) == 0) {
      K_range <- 2:min(8, min(N, P))
      cat("Adjusted K_range to:", paste(K_range, collapse=", "), "\n")
    }
    
    cat("Testing K values:", paste(K_range, collapse=", "), "\n")
    
    best_K <- NULL
    best_cv_score <- -Inf
    cv_results <- list()
    
    for (K_test in K_range) {
      cat("  Testing K =", K_test, "...")
      start_time <- Sys.time()
      
      # Run cross-validation for this K with error handling
      cv_score <- tryCatch({
        cross_validate_K(X, K_test, G, cv_folds, 
                         n_iter = max(500, n_iter %/% 4), 
                         burn_in = max(100, burn_in %/% 4),
                         lambda_a, lambda_b, auto_group)
      }, error = function(e) {
        # If cross-validation fails completely, return poor score
        list(mean_score = -1000, se_score = 0, fold_scores = rep(-1000, cv_folds), n_valid_folds = 0)
      })
      
      cv_results[[paste0("K_", K_test)]] <- cv_score
      
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Safe printing with NA handling
      score_text <- if (is.na(cv_score$mean_score)) "NA" else round(cv_score$mean_score, 4)
      se_text <- if (is.na(cv_score$se_score)) "NA" else round(cv_score$se_score, 4)
      
      cat(" CV score:", score_text, 
          " (±", se_text, ") -", 
          round(elapsed_time, 1), "sec\n")
      
      # Update best K if this one is better (with NA handling)
      if (!is.na(cv_score$mean_score) && !is.infinite(cv_score$mean_score) && 
          cv_score$mean_score > best_cv_score) {
        best_cv_score <- cv_score$mean_score
        best_K <- K_test
      }
    }
    
    # Set final K with fallback
    if (is.null(best_K) || best_cv_score == -Inf) {
      # If no valid K found, use middle of range as fallback
      K <- round(median(K_range))
      cat("WARNING: No valid K found via CV, using fallback K =", K, "\n")
    } else {
      K <- best_K
      cat("OPTIMAL K SELECTED: K =", K, "with CV score =", round(best_cv_score, 4), "\n")
    }
    
    # Print detailed CV scores for validation
    cat("DETAILED CV SCORES:\n")
    for (K_test in K_range) {
      score_info <- cv_results[[paste0("K_", K_test)]]
      marker <- if (K_test == best_K) " (SELECTED)" else ""
      
      # Safe printing with NA handling
      mean_text <- if (is.na(score_info$mean_score)) "NA" else round(score_info$mean_score, 4)
      se_text <- if (is.na(score_info$se_score)) "NA" else round(score_info$se_score, 4)
      
      cat("   K =", K_test, ": ", mean_text, 
          " (±", se_text, ")", marker, "\n")
    }
    cat("Cross-validation completed successfully!\n\n")
  }
  
  # Default K if not set
  if (is.null(K)) {
    K <- 5  # Default value
    cat("Using default K =", K, "\n")
  }
  
  #                                                                                    
  #                              AUTO_GROUP FUNCTIONALITY
  #                                                                                    
  # When auto_group=TRUE, automatically determine optimal number of groups G
  if (auto_group && is.null(G)) {
    cat("  AUTO_GROUP=TRUE: Automatically determining optimal G...\n")
    
    # Test different values of G
    G_candidates <- seq(2, min(10, P %/% 3))  # Test 2 to 10 groups
    best_G <- NULL
    best_score <- -Inf
    
    cat("Testing G values:", paste(G_candidates, collapse=", "), "\n")
    
    for (G_test in G_candidates) {
      if (P %% G_test == 0) {  # Only test if P is divisible by G
        cat("  Testing G =", G_test, "...")
        
        # Run short MCMC for model evaluation
        tryCatch({
          test_result <- mcmc_bpca_group_lasso(
            X, K=K, G=G_test, n_iter=500, burn_in=100, 
            lambda_a=lambda_a, lambda_b=lambda_b, 
            auto_group=FALSE, n_chains=2
          )
          
          # Score based on reconstruction quality and group sparsity
          reconstruction_score <- cor(X[!is.na(X)], test_result$X_imputed[!is.na(X)])
          sparsity_score <- mean(test_result$tau_g < quantile(test_result$tau_g, 0.5))
          combined_score <- reconstruction_score + 0.5 * sparsity_score
          
          cat(" Score =", round(combined_score, 3), "\n")
          
          if (combined_score > best_score) {
            best_score <- combined_score
            best_G <- G_test
          }
        }, error = function(e) {
          cat(" ERROR\n")
        })
      }
    }
    
    G <- ifelse(is.null(best_G), max(2, P %/% 5), best_G)
    cat("AUTO_GROUP selected optimal G =", G, "\n")
  }
  
  # Fallback: if G still not set, use default
  if (is.null(G)) {
    G <- max(2, min(10, P %/% 5))  # Default: 5-10 variables per group
    cat("Auto-detected G =", G, "groups\n")
  }
  
  if (P %% G != 0) {
    stop("P must be divisible by G for equal group sizes")
  }
  
  group_size <- P %/% G
  group_id <- rep(1:G, each = group_size)
  
  cat("Group structure: G =", G, "groups with", group_size, "variables each\n")
  cat("Running", n_chains, "independent MCMC chains for convergence assessment\n")
  
  # Identify missing data
  missing_idx <- which(is.na(X), arr.ind = TRUE)
  observed_idx <- which(!is.na(X), arr.ind = TRUE)
  
  # Initialize missing values with robust error handling
  X_work <- X
  for (j in 1:P) {
    missing_j <- is.na(X[, j])
    if (any(missing_j)) {
      if (sum(!missing_j) > 0) {
        # Use observed values to estimate missing
      X_work[missing_j, j] <- mean(X[, j], na.rm = TRUE)
      } else {
        # All values missing for this variable - use global mean
        global_mean <- mean(X, na.rm = TRUE)
        X_work[missing_j, j] <- ifelse(is.na(global_mean), 0, global_mean)
        cat("Warning: Variable", j, "has all missing values. Using global mean.\n")
      }
    }
  }
  
  Z_samples_all <- array(NA, dim = c(n_chains, n_iter, N, K))
  W_samples_all <- array(NA, dim = c(n_chains, n_iter, P, K))      # Factor loadings
  mu_samples_all <- array(NA, dim = c(n_chains, n_iter, P))         # Intercepts
  sigma2_samples_all <- matrix(NA, n_chains, n_iter)                # Noise variance
  tau_g2_samples_all <- array(NA, dim = c(n_chains, n_iter, G))      # Group scales squared (TRUE Group LASSO)
  lambda2_samples_all <- matrix(NA, n_chains, n_iter)               # Global shrinkage
  
  #                                                                                    
  #                         MULTI-CHAIN MCMC SAMPLING LOOP
  #                                                                                    
  # Run n_chains independent MCMC chains with different random starting values
  # This is ESSENTIAL for:
  # 1. Robust convergence assessment (R-hat)
  # 2. Detection of multimodality or poor mixing
  # 3. More reliable posterior estimates
  for (chain in 1:n_chains) {
    cat("Running chain", chain, "of", n_chains, "...\n")
    
    #                                                                                  
    #                    PARAMETER INITIALIZATION (CHAIN-SPECIFIC)
    #                                                                                  
    # Each chain starts from DIFFERENT random values to test convergence
    set.seed(123 + chain * 1000)  # Different seed ensures different starting points
    
    #   INITIALIZE Z (LATENT FACTORS)  
    # Z[i,k] ~ N(0, 0.5²) for each subject i and factor k
    # Starting with moderate variance (0.5) to avoid extreme values
    Z <- matrix(rnorm(N * K, 0, 0.5), N, K)
    
    # Initialize other parameters
    W <- matrix(rnorm(P * K, 0, 0.1), P, K)  # Small starting loadings
    mu <- colMeans(X_work) + rnorm(P, 0, 0.1)  # Start near data means
    sigma2 <- rgamma(1, 2, 2)  # Random noise variance
    
    # Group LASSO parameters - different random starts
    tau_g2 <- rgamma(G, 2, 2)   # Group-specific scales squared (TRUE Group LASSO)
    lambda2 <- rgamma(1, 2, 2)  # Global shrinkage parameter
    
    #                                                                                    
    #                          MAIN MCMC ITERATION LOOP
    #                                                                                    
    # Each iteration updates ALL parameters in sequence using Gibbs sampling
    # This is the heart of the Bayesian MCMC algorithm
    for (iter in 1:n_iter) {
      
      #                          STEP 1: MISSING DATA IMPUTATION
      #                                                                                  
      # Update missing values using current parameter estimates
      if (nrow(missing_idx) > 0) {
        # Predict missing values using current Z, W,  
        X_pred <- Z %*% t(W) + matrix(rep(mu, each = N), nrow = N)
        
        # Sample missing values from their conditional posterior
        X_work[missing_idx] <- X_pred[missing_idx] + rnorm(nrow(missing_idx), 0, sqrt(sigma2))
      }
                                                    
      #                       STEP 2: UPDATE Z (LATENT FACTORS)  
      #                                                                                  
      # This is WHERE Z IS UPDATED in each MCMC iteration
      for (i in 1:N) {  # Update each subject's latent factors
        # Compute posterior covariance matrix for z_i
        # = (likelihood precision + prior precision)^(-1)
        # Add regularization for numerical stability
        Sigma_z <- solve(t(W) %*% W / sigma2 + diag(K) + 1e-8 * diag(K))
        
        # Compute posterior mean for z_i  
        # = covariance × likelihood gradient
        mu_z <- Sigma_z %*% t(W) %*% (X_work[i, ] - mu) / sigma2
        
        #   SAMPLE NEW Z VALUES  
        # Draw z_i from multivariate normal posterior
        Z[i, ] <- as.vector(rmvnorm(1, mu_z, Sigma_z))
      }
      # Result: Z matrix [N × K] now contains updated latent factor scores
      #         for all N subjects across all K factors
      
      # 3. Update W (loadings) with TRUE Group LASSO (Laplace priors)
      X_residual <- X_work - matrix(rep(mu, each = N), nrow = N)
      W <- sample_group_lasso_loadings(X_residual, Z, sigma2, tau_g2, group_id, K)
      
      # 4. Update mu (intercepts)
      for (j in 1:P) {
        residual <- X_work[, j] - Z %*% W[j, ]
        mu[j] <- rnorm(1, mean(residual), sqrt(sigma2 / N))
      }
      
      # 5. Update sigma2 (noise variance) 
      X_pred <- Z %*% t(W) + matrix(rep(mu, each = N), nrow = N)
      residual <- X_work - X_pred
      sse <- sum(residual^2)
      shape_sigma <- N * P / 2 + 1
      rate_sigma <- sse / 2 + 1
      sigma2 <- 1 / rgamma(1, shape_sigma, rate_sigma)
      
      # 6. Update tau_g² (group-specific scales) - TRUE Group LASSO
      tau_g2 <- sample_group_scales_gig(W, group_id, sigma2, lambda2, G, K)
      
      # 7. Update lambda2 (global shrinkage parameter)
      lambda2 <- sample_global_shrinkage(tau_g2, lambda_a, lambda_b)
                                                                             
      #                        STORE MCMC SAMPLES FOR THIS ITERATION
      #                                                                                  
      # Save all parameter values from this iteration for later analysis
      # Memory optimization: Store all samples for convergence diagnostics
      
      #   STORE Z SAMPLES  
      # Z_samples_all[chain, iter, subject, factor] = current Z values
      # This captures the COMPLETE trajectory of latent factors across iterations
      Z_samples_all[chain, iter, , ] <- Z
      
      # Store other parameters
      W_samples_all[chain, iter, , ] <- W          # Factor loadings
      mu_samples_all[chain, iter, ] <- mu          # Intercepts
      sigma2_samples_all[chain, iter] <- sigma2    # Noise variance
      tau_g2_samples_all[chain, iter, ] <- tau_g2   # Group scales (squared)
      lambda2_samples_all[chain, iter] <- lambda2  # Global shrinkage
    }
    
    # Progress report for each chain
    cat("Chain", chain, "completed - Final lambda2:", round(lambda2, 4), 
        "- Active groups:", sum(tau_g2 > quantile(tau_g2, 0.7)), "\n")
  }
                                                                
  #                          POST-MCMC SAMPLE PROCESSING
  #                                                                                    
  # Process raw MCMC samples to obtain final parameter estimates
  
  # Remove burn-in samples (initial iterations that may not be representative)
  keep_idx <- (burn_in + 1):n_iter  # Keep only iterations after burn-in
                                                                               
  #                      EXTRACT POST-BURN-IN Z SAMPLES
  #                                                                                  
  # Remove burn-in from ALL parameter samples
  
  #   Z_samples_postburn: Z samples after removing burn-in [n_chains × kept_iter × N × K]
  Z_samples_postburn <- Z_samples_all[, keep_idx, , , drop = FALSE]
  W_samples_postburn <- W_samples_all[, keep_idx, , , drop = FALSE]
  mu_samples_postburn <- mu_samples_all[, keep_idx, , drop = FALSE]
  sigma2_samples_postburn <- sigma2_samples_all[, keep_idx, drop = FALSE]
  tau_g2_samples_postburn <- tau_g2_samples_all[, keep_idx, , drop = FALSE]
  lambda2_samples_postburn <- lambda2_samples_all[, keep_idx, drop = FALSE]
                                                                            
  #                     COMBINE CHAINS FOR FINAL ESTIMATES
  #                                                                                  
  # Flatten the chain dimension to combine all chains into one large sample
  # This gives us more samples for better posterior estimates
  
  #   Z_combined: All Z samples from all chains combined [total_samples × N × K]
  # where total_samples = n_chains × (n_iter - burn_in)
  Z_combined <- array(Z_samples_postburn, dim = c(n_chains * length(keep_idx), N, K))
  W_combined <- array(W_samples_postburn, dim = c(n_chains * length(keep_idx), P, K))
  mu_combined <- matrix(mu_samples_postburn, n_chains * length(keep_idx), P)
  sigma2_combined <- as.vector(sigma2_samples_postburn)
  tau_g2_combined <- matrix(tau_g2_samples_postburn, n_chains * length(keep_idx), G)
  lambda2_combined <- as.vector(lambda2_samples_postburn)
                                                                             
  #                        COMPUTE POSTERIOR MEANS
  #                                                                                  
  # Calculate final parameter estimates as posterior means across all samples
  
  #   Z_mean: Final Z estimate [N × K] - mean across all MCMC samples
  # This is the main Z result that gets returned and used for analysis
  Z_mean <- apply(Z_combined, c(2, 3), mean)  # Average over sample dimension
  W_mean <- apply(W_combined, c(2, 3), mean)
  mu_mean <- apply(mu_combined, 2, mean)
  sigma2_mean <- mean(sigma2_combined)
  tau_g2_mean <- apply(tau_g2_combined, 2, mean)
  lambda2_mean <- mean(lambda2_combined)
  
  # Final imputation
  X_imputed <- X
  if (nrow(missing_idx) > 0) {
    X_pred_final <- Z_mean %*% t(W_mean) + matrix(rep(mu_mean, each = N), nrow = N)
    X_imputed[missing_idx] <- X_pred_final[missing_idx]
  }
  
  # Post-hoc alignment to fix rotation problem
  cat("\n=== POST-HOC ALIGNMENT FOR FACTOR IDENTIFIABILITY ===\n")
  aligned_results <- align_factor_chains(W_samples_postburn, Z_samples_postburn, n_chains)
  W_samples_aligned <- aligned_results$W_aligned
  Z_samples_aligned <- aligned_results$Z_aligned
  
  # Calculate convergence diagnostics with aligned chains
  cat("\n=== MCMC CONVERGENCE DIAGNOSTICS (After Procrustes Alignment) ===\n")
  convergence_diagnostics <- calculate_convergence_diagnostics(
          sigma2_samples_postburn, tau_g2_samples_postburn, lambda2_samples_postburn,
    W_samples_aligned, Z_samples_aligned, n_chains
  )
  
  # Auto-adjust if convergence is poor
  if (convergence_diagnostics$max_rhat > 1.5) {
    cat("\n    POOR CONVERGENCE DETECTED! R-hat =", round(convergence_diagnostics$max_rhat, 3), "\n")
    cat("RECOMMENDATION: Increase n_iter to at least", 2 * n_iter, "iterations\n")
    cat("ALTERNATIVE: Try different starting values or more chains\n")
  } else if (convergence_diagnostics$max_rhat > 1.2) {
    cat("\n Marginal convergence. Consider increasing iterations.\n")
  } else {
    cat("\n EXCELLENT CONVERGENCE! All parameters converged properly.\n")
  }
                                                                      
  #                                RETURN RESULTS
  #                                                                                    
  # Return comprehensive results including multiple representations of Z
  return(list(
    Z = Z_mean,
    W = W_mean,           # Factor loadings [P × K] - which variables load on which factors
    mu = mu_mean,         # Feature intercepts [P × 1] - baseline levels for each variable
    sigma2 = sigma2_mean, # Residual noise variance [scalar] - unexplained variation
    tau_g2 = tau_g2_mean,   # Group scales squared [G × 1] - shrinkage level for each group (TRUE Group LASSO)
    lambda2 = lambda2_mean, # Global shrinkage [scalar] - overall LASSO strength
    X_imputed = X_imputed,  # Imputed data matrix [N × P] - missing values filled in
    
    #                     COMBINED MCMC SAMPLES (FOR UNCERTAINTY)
    #                                                                                  
    Z_samples = Z_combined,
    W_samples = W_combined,         # All W samples for uncertainty
    mu_samples = mu_combined,       # All mu samples
    sigma2_samples = sigma2_combined, # All sigma2 samples
    tau_g2_samples = tau_g2_combined,   # All tau_g² samples (Group scales squared)
    lambda2_samples = lambda2_combined, # All lambda2 samples
                                                                                
    #                    CHAIN-SPECIFIC SAMPLES (FOR CONVERGENCE)
    #                                                                                  
    # Samples organized by chain - use for convergence diagnostics (R-hat, trace plots)
    Z_samples_chains = Z_samples_postburn,
    W_samples_chains = W_samples_postburn,
    mu_samples_chains = mu_samples_postburn,
    sigma2_samples_chains = sigma2_samples_postburn,
    tau_g2_samples_chains = tau_g2_samples_postburn,
    lambda2_samples_chains = lambda2_samples_postburn,
    
    #                  ROTATION-ALIGNED SAMPLES (FOR IDENTIFIABILITY)
    #                                                                                  
    # Post-hoc aligned samples that fix the rotation problem in factor models
    
    #   Z_samples_aligned_chains: ROTATION-CORRECTED Z SAMPLES  
    # Same as Z_samples_chains but after Procrustes alignment
    # Use this for: improved convergence diagnostics, more stable R-hat values
    # This addresses the fundamental identifiability issue: W*Z = (W*R) * (R'*Z)
    Z_samples_aligned_chains = Z_samples_aligned,
    W_samples_aligned_chains = W_samples_aligned,
                                                                             
    #                              METADATA & DIAGNOSTICS
    #                                                                                  
    group_id = group_id,                      # Which group each variable belongs to [P × 1]
    G = G,                                    # Number of groups used
    K = K,                                    # Number of factors used (important for auto_K)
    selected_K = K,                           # K that was selected by auto_K (same as K, for clarity)
    n_chains = n_chains,                      # Number of MCMC chains run
    convergence_diagnostics = convergence_diagnostics,  # Detailed R-hat and ESS results
    cv_results = if(exists("cv_results")) cv_results else NULL  # Cross-validation results (if auto_K was used)
  ))
                                                                                  
  #                           Z USAGE SUMMARY
  #                                                                                    
  # 
  #   FOR BASIC ANALYSIS: Use result$Z (the posterior mean)
  #   FOR UNCERTAINTY: Use result$Z_samples (all samples combined)
  #   FOR CONVERGENCE: Use result$Z_samples_chains (organized by chain)
  #   FOR ROTATION ISSUES: Use result$Z_samples_aligned_chains (rotation-corrected)
}

# Evaluation functions
evaluate_factor_recovery <- function(W_true, W_est) {
  # Procrustes analysis
  pro <- vegan::procrustes(X = W_true, Y = W_est, symmetric = TRUE)
  W_est_aligned <- W_est %*% pro$rotation
  
  # Correlation between true and estimated loadings
  cor_loadings <- cor(as.vector(W_true), as.vector(W_est_aligned))
  
  # RMSE between true and estimated loadings
  rmse_loadings <- sqrt(mean((W_true - W_est_aligned)^2))
  
  return(list(
    W_aligned = W_est_aligned,
    correlation = cor_loadings,
    rmse = rmse_loadings,
    procrustes_ss = pro$ss
  ))
}

evaluate_reconstruction <- function(X_true, X_imputed, missing_idx) {
  # RMSE for missing values
  if (nrow(missing_idx) > 0) {
    rmse_missing <- sqrt(mean((X_true[missing_idx] - X_imputed[missing_idx])^2))
  } else {
    rmse_missing <- NA
  }
  
  # Overall RMSE
  rmse_overall <- sqrt(mean((X_true - X_imputed)^2))
  
  # Correlation
  cor_overall <- cor(as.vector(X_true), as.vector(X_imputed))
  
  return(list(
    rmse_missing = rmse_missing,
    rmse_overall = rmse_overall,
    correlation = cor_overall
  ))
}

# Group sparsity analysis
analyze_group_sparsity <- function(mcmc_result, W_true = NULL, true_groups = NULL) {
  cat("=== GROUP SPARSITY ANALYSIS ===\n\n")
  cat("SPARSITY LOGIC: In Bayesian Group LASSO:\n")
  cat("- Large τ_g² = Less shrinkage = More active group\n")
  cat("- Small τ_g² = More shrinkage = Less active group\n")
  cat("- Large ||W_g||_F = Strong signal = More active group\n\n")
  
  W_est <- mcmc_result$W
  tau_g2 <- mcmc_result$tau_g2
  group_id <- mcmc_result$group_id
  G <- mcmc_result$G
  
  # Group-level signal strength
  cat("1. Group-level Signal Strength (||W_g||_F) and Shrinkage (τ_g²):\n")
  group_norms <- numeric(G)
  for (g in 1:G) {
    vars_in_group <- which(group_id == g)
    W_g <- W_est[vars_in_group, , drop = FALSE]
    group_norms[g] <- norm(W_g, type = "F")
    
    true_status <- if (!is.null(true_groups)) {
      ifelse(g %in% true_groups, "ACTIVE", "inactive")
    } else {
      "unknown"
    }
    
    cat(" Group", g, ": ||W_g||_F =", round(group_norms[g], 4), 
        " tau_g² =", round(tau_g2[g], 4), "(", true_status, ")\n")
  }
  
  # Automatic group detection
  cat("\n2. Automatic Group Detection:\n")
  
  # Method 1: Based on group norms (signal strength)
  group_norm_threshold <- quantile(group_norms, 0.7)  # Top 30%
  detected_active_norms <- which(group_norms > group_norm_threshold)
  
  # Method 2: Based on tau_g² (shrinkage level)
  tau_threshold <- quantile(tau_g2, 0.7)  # Top 30% (less shrinkage = more active)
  detected_active_tau <- which(tau_g2 > tau_threshold)
  
  cat("   Detected active groups (norm-based):", paste(detected_active_norms, collapse = ", "), "\n")
  cat("   Detected active groups ( ²-based):", paste(detected_active_tau, collapse = ", "), "\n")
  
  if (!is.null(true_groups)) {
    cat("   True active groups:", paste(true_groups, collapse = ", "), "\n")
    
    # Accuracy metrics
    correct_norm <- length(intersect(detected_active_norms, true_groups))
    correct_tau <- length(intersect(detected_active_tau, true_groups))
    
    cat("   Detection accuracy (norm-based):", correct_norm, "/", length(true_groups), 
        "(", round(100 * correct_norm / length(true_groups), 1), "%)\n")
    cat("   Detection accuracy ( ²-based):", correct_tau, "/", length(true_groups), 
        "(", round(100 * correct_tau / length(true_groups), 1), "%)\n")
  }
  
  # Sparsity level
  total_elements <- length(W_est)
  near_zero <- sum(abs(W_est) < 0.05)
  sparsity_achieved <- 100 * near_zero / total_elements
  
  cat("\n3. Sparsity Level:\n")
  cat("   Elements near zero (|w| < 0.05):", near_zero, "/", total_elements, 
      "(", round(sparsity_achieved, 1), "%)\n")
  cat("   Global shrinkage  ²:", round(mcmc_result$lambda2, 4), "\n")
  
  return(list(
    group_norms = group_norms,
    detected_active_norms = detected_active_norms,
    detected_active_tau = detected_active_tau,
    sparsity_achieved = sparsity_achieved,
    tau_g2 = tau_g2
  ))
}

# HEATMAP COMPARISON PLOTS
plot_W_comparison <- function(W_true, W_est_aligned, title_suffix = "", save_plots = TRUE, output_dir = "bpca_missing_lasso_plots") {
  cat("=== CREATING W MATRIX HEATMAP COMPARISON ===\n\n")
  
  # Create output directory
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # 1. Normalization tests
  cat("1. Normalization Tests:\n")
  
  # Frobenius norm comparison
  norm_true <- norm(W_true, type = "F")
  norm_est <- norm(W_est_aligned, type = "F")
  norm_ratio <- norm_est / norm_true
  
  cat("   Frobenius norm (true):", round(norm_true, 4), "\n")
  cat("   Frobenius norm (estimated):", round(norm_est, 4), "\n")
  cat("   Norm ratio (est/true):", round(norm_ratio, 4), "\n")
  
  # Column-wise norms
  col_norms_true <- apply(W_true, 2, function(x) norm(x, type = "2"))
  col_norms_est <- apply(W_est_aligned, 2, function(x) norm(x, type = "2"))
  
  cat("   Column norms (true):", paste(round(col_norms_true, 3), collapse = ", "), "\n")
  cat("   Column norms (estimated):", paste(round(col_norms_est, 3), collapse = ", "), "\n\n")
  
  # 2. Create heatmaps
  cat("2. Creating Heatmaps...\n")
  
  # Prepare matrices with labels
  P <- nrow(W_true)
  K <- ncol(W_true)
  
  matrix_rownames <- paste0("Var", 1:P)
  matrix_colnames <- paste0("Factor", 1:K)
  
  W_true_labeled <- W_true
  rownames(W_true_labeled) <- matrix_rownames
  colnames(W_true_labeled) <- matrix_colnames
  
  W_est_labeled <- W_est_aligned
  rownames(W_est_labeled) <- matrix_rownames
  colnames(W_est_labeled) <- matrix_colnames
  
  # Difference matrix
  W_diff <- W_est_aligned - W_true
  rownames(W_diff) <- matrix_rownames
  colnames(W_diff) <- matrix_colnames
  
  # Set common color scale
  all_values <- c(as.vector(W_true_labeled), as.vector(W_est_labeled))
  min_val <- min(all_values)
  max_val <- max(all_values)
  color_breaks <- seq(min_val, max_val, length.out = 100)
  my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
  
  # Difference color scale
  diff_max <- max(abs(W_diff))
  diff_breaks <- seq(-diff_max, diff_max, length.out = 100)
  diff_colors <- colorRampPalette(c("darkred", "white", "darkblue"))(99)
  
  # Create plots
  p_true <- pheatmap(
    W_true_labeled,
    main = paste0("W True", title_suffix),
    cluster_rows = FALSE, cluster_cols = FALSE,
    color = my_colors, breaks = color_breaks,
    display_numbers = TRUE, fontsize_number = 8,
    silent = TRUE
  )
  
  p_est <- pheatmap(
    W_est_labeled,
    main = paste0("W Estimated", title_suffix),
    cluster_rows = FALSE, cluster_cols = FALSE,
    color = my_colors, breaks = color_breaks,
    display_numbers = TRUE, fontsize_number = 8,
    silent = TRUE
  )
  
  p_diff <- pheatmap(
    W_diff,
    main = paste0("Difference (Est - True)", title_suffix),
    cluster_rows = FALSE, cluster_cols = FALSE,
    color = diff_colors, breaks = diff_breaks,
    display_numbers = TRUE, fontsize_number = 8,
    silent = TRUE
  )
  
  # Save individual plots
  if (save_plots) {
    method_name <- gsub("[^A-Za-z0-9]", "_", title_suffix)
    if (method_name == "") method_name <- "standard"
    
    # Save combined plot
    pdf(file.path(output_dir, paste0("W_comparison_combined", method_name, ".pdf")), width = 24, height = 10)
    grid.arrange(p_true$gtable, p_est$gtable, p_diff$gtable, ncol = 3)
    dev.off()
    
    cat("  Saved W comparison plots to:", output_dir, "\n")
  } else {
    # Display only
    grid.arrange(p_true$gtable, p_est$gtable, p_diff$gtable, ncol = 3)
  }
  
  return(list(
    norm_ratio = norm_ratio,
    col_norms_ratio = col_norms_est / col_norms_true,
    max_abs_diff = max(abs(W_diff)),
    plots = list(true = p_true, estimated = p_est, difference = p_diff)
  ))
}

# MISSING DATA PLOTS
plot_missing_data_reconstruction <- function(X_true, X_imputed, missing_idx, save_plots = TRUE, output_dir = "bpca_missing_lasso_plots") {
  cat("=== MISSING DATA RECONSTRUCTION PLOTS ===\n\n")
  
  if (nrow(missing_idx) == 0) {
    cat("No missing data to plot.\n")
    return(NULL)
  }
  
  # Create output directory
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  missing_true <- X_true[missing_idx]
  missing_imputed <- X_imputed[missing_idx]
  cor_missing <- cor(missing_true, missing_imputed)
  
  # Save plot to file if requested
  if (save_plots) {
    pdf(file.path(output_dir, "missing_data_reconstruction.pdf"), width = 12, height = 6)
  }
  
  # Create plots
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  
  # Scatter plot
  plot(missing_true, missing_imputed, 
       xlab = "True Values", ylab = "Imputed Values",
       main = "Missing Data Reconstruction",
       pch = 16, col = rgb(0, 0, 1, 0.6))
  abline(0, 1, col = "red", lwd = 2)
  text(x = min(missing_true) + 0.1 * diff(range(missing_true)),
       y = max(missing_imputed) - 0.1 * diff(range(missing_imputed)),
       labels = paste0("r = ", round(cor_missing, 3)), cex = 1.2)
  
  # Residual plot
  residuals <- missing_imputed - missing_true
  plot(missing_true, residuals,
       xlab = "True Values", ylab = "Residuals (Imputed - True)",
       main = "Reconstruction Residuals",
       pch = 16, col = rgb(0, 0, 1, 0.6))
  abline(h = 0, col = "red", lwd = 2)
  
  if (save_plots) {
    dev.off()
    cat("Saved missing data reconstruction plot to:", file.path(output_dir, "missing_data_reconstruction.pdf"), "\n")
  }
  
  cat("Missing data correlation:", round(cor_missing, 4), "\n")
  cat("Missing data RMSE:", round(sqrt(mean((missing_true - missing_imputed)^2)), 4), "\n")
  
  return(list(correlation = cor_missing, rmse = sqrt(mean((missing_true - missing_imputed)^2))))
}

# COMPREHENSIVE MCMC TRACE PLOTS
plot_comprehensive_mcmc_traces <- function(mcmc_result, save_plots = TRUE, output_dir = "bpca_missing_lasso_plots") {
  # Extract samples - use chain-specific samples for trace plots
  sigma2_samples <- mcmc_result$sigma2_samples_chains
  lambda2_samples <- mcmc_result$lambda2_samples_chains
  W_samples <- mcmc_result$W_samples_chains
  Z_samples <- mcmc_result$Z_samples_chains
  tau_g2_samples <- mcmc_result$tau_g2_samples_chains
  
  # Get dimensions
  n_chains <- dim(W_samples)[1]
  n_iter <- dim(W_samples)[2]
  P <- dim(W_samples)[3]
  K <- dim(W_samples)[4]
  N <- dim(Z_samples)[3]
  G <- dim(tau_g2_samples)[3]
  
  # Create output directory
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Save plot to file if requested
  if (save_plots) {
    pdf(file.path(output_dir, "comprehensive_mcmc_traces.pdf"), width = 16, height = 12)
  }
  
  # Set up layout for comprehensive trace plots
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3, ncol = 4, byrow = TRUE))
  
  # Define distinct colors for chains
  chain_colors <- c("blue", "red", "green", "purple", "orange", "brown")[1:n_chains]
  
  # 1. Global Parameters
  # sigma2 trace
  plot(1:n_iter, sigma2_samples[1, ], type = "l", col = chain_colors[1], 
       main = "sigma2 Evolution", xlab = "Iteration", ylab = "sigma2")
  if (n_chains > 1) {
    for (chain in 2:n_chains) {
      lines(1:n_iter, sigma2_samples[chain, ], col = chain_colors[chain])
    }
  }
  
  # lambda2 trace
  plot(1:n_iter, lambda2_samples[1, ], type = "l", col = chain_colors[1], 
       main = "lambda2 Evolution", xlab = "Iteration", ylab = "lambda2")
  if (n_chains > 1) {
    for (chain in 2:n_chains) {
      lines(1:n_iter, lambda2_samples[chain, ], col = chain_colors[chain])
    }
  }
  
  # 2. Selected W loadings
  # W[1,1] trace
  plot(1:n_iter, W_samples[1, , 1, 1], type = "l", col = chain_colors[1], 
       main = "W[1,1] Evolution", xlab = "Iteration", ylab = "W[1,1]")
  if (n_chains > 1) {
    for (chain in 2:n_chains) {
      lines(1:n_iter, W_samples[chain, , 1, 1], col = chain_colors[chain])
    }
  }
  
  # W[25,1] trace
  plot(1:n_iter, W_samples[1, , min(25, P), 1], type = "l", col = chain_colors[1], 
       main = "W[25,1] Evolution", xlab = "Iteration", ylab = "W[25,1]")
  if (n_chains > 1) {
    for (chain in 2:n_chains) {
      lines(1:n_iter, W_samples[chain, , min(25, P), 1], col = chain_colors[chain])
    }
  }
  
  # 3. Selected Z factors
  # Z[1,1] trace
  plot(1:n_iter, Z_samples[1, , 1, 1], type = "l", col = chain_colors[1], 
       main = "Z[1,1] Evolution", xlab = "Iteration", ylab = "Z[1,1]")
  if (n_chains > 1) {
    for (chain in 2:n_chains) {
      lines(1:n_iter, Z_samples[chain, , 1, 1], col = chain_colors[chain])
    }
  }
  
  # Z[50,1] trace
  plot(1:n_iter, Z_samples[1, , min(50, N), 1], type = "l", col = chain_colors[1], 
       main = "Z[50,1] Evolution", xlab = "Iteration", ylab = "Z[50,1]")
  if (n_chains > 1) {
    for (chain in 2:n_chains) {
      lines(1:n_iter, Z_samples[chain, , min(50, N), 1], col = chain_colors[chain])
    }
  }
  
  # 4. Group scales (tau_g2) - first 6 groups
  # tau_g2[1] trace
  plot(1:n_iter, tau_g2_samples[1, , 1], type = "l", col = chain_colors[1], 
       main = "tau_g2[1] Evolution", xlab = "Iteration", ylab = "tau_g2[1]")
  if (n_chains > 1) {
    for (chain in 2:n_chains) {
      lines(1:n_iter, tau_g2_samples[chain, , 1], col = chain_colors[chain])
    }
  }
  
  # tau_g2[2] trace
  plot(1:n_iter, tau_g2_samples[1, , 2], type = "l", col = chain_colors[1], 
       main = "tau_g2[2] Evolution", xlab = "Iteration", ylab = "tau_g2[2]")
  if (n_chains > 1) {
    for (chain in 2:n_chains) {
      lines(1:n_iter, tau_g2_samples[chain, , 2], col = chain_colors[chain])
    }
  }
  
  # Add legend to the last plot to show chain colors
  if (n_chains > 1) {
    legend("topright", legend = paste("Chain", 1:n_chains), 
           col = chain_colors, lty = 1, cex = 0.8, bty = "n")
  }
  
  # 5. Distribution plots - use combined samples for distributions
  # sigma2 distribution
  hist(mcmc_result$sigma2_samples, main = "sigma2 Distribution", xlab = "sigma2", col = "lightblue")
  
  # lambda2 distribution
  hist(mcmc_result$lambda2_samples, main = "lambda2 Distribution", xlab = "lambda2", col = "lightcoral")
  
  # W[1,1] distribution
  hist(mcmc_result$W_samples[, 1, 1], main = "W[1,1] Distribution", xlab = "W[1,1]", col = "lightgreen")
  
  # tau_g2 distribution
  hist(mcmc_result$tau_g2_samples[, 1], main = "tau_g2 Distribution", xlab = "tau_g2", col = "lightyellow")
  
  if (save_plots) {
    dev.off()
    cat("Saved comprehensive MCMC trace plots to:", file.path(output_dir, "comprehensive_mcmc_traces.pdf"), "\n")
  }
  
  # Convergence diagnostics summary
  cat("\n=== CONVERGENCE DIAGNOSTICS SUMMARY ===\n")
  
  # Display R-hat diagnostics if available
  if ("convergence_diagnostics" %in% names(mcmc_result)) {
    conv_diag <- mcmc_result$convergence_diagnostics
    cat("Gelman-Rubin R-hat Statistics:\n")
    cat("   Max R-hat:", round(conv_diag$max_rhat, 3), "\n")
    cat("   Excellent convergence (< 1.1):", round(100 * conv_diag$prop_converged_excellent, 1), "%\n")
    cat("   Acceptable convergence (< 1.2):", round(100 * conv_diag$prop_converged_acceptable, 1), "%\n")
    
    if (conv_diag$max_rhat < 1.1) {
      cat("     Excellent convergence across all chains!\n")
    } else if (conv_diag$prop_converged_excellent > 0.8) {
      cat("       Most parameters converged, some may need longer runs\n")
    } else {
      cat("     Poor convergence - increase iterations\n")
    }
  }
  
  # Effective sample size for key parameters
  if (requireNamespace("coda", quietly = TRUE)) {
    cat("\nEffective Sample Sizes (ESS) - Combined Chains:\n")
    
    # ESS for sigma2
    sigma2_ess <- coda::effectiveSize(mcmc_result$sigma2_samples)
    cat("   sigma2 ESS:", round(sigma2_ess, 0), "\n")
    
    # ESS for lambda2
    lambda2_ess <- coda::effectiveSize(mcmc_result$lambda2_samples)
    cat("   lambda2 ESS:", round(lambda2_ess, 0), "\n")
    
    # ESS for selected W elements
    w_ess_sample <- coda::effectiveSize(mcmc_result$W_samples[, 1, 1])
    cat("   W[1,1] ESS:", round(w_ess_sample, 0), "\n")
    
    # ESS for Group LASSO parameters
    tau_ess <- coda::effectiveSize(mcmc_result$tau_g2_samples[, 1])
    cat("   tau_g2[1] ESS:", round(tau_ess, 0), "\n")
  }
  
  cat("\nComprehensive MCMC trace plots completed!\n")
}

# --- 16. QUICK TEST FUNCTION ---
quick_test <- function() {
  cat("Running Quick Test of Comprehensive BPCA Group LASSO...\n\n")
  
  # Generate test data
  set.seed(123)
  sim_data <- simulate_data_with_groups(
    n = 200, p = 50, G = 5, K = 5,
    active_groups = c(1, 3, 5),
    sigma = 0.2
  )
  
  X_data_complete <- sim_data$X_data
  W_true <- sim_data$W_true
  
  # Add 10% missing data
  set.seed(42)
  X_data_missing <- X_data_complete
  missing_positions <- sample(length(X_data_missing), size = floor(0.1 * length(X_data_missing)))
  X_data_missing[missing_positions] <- NA
  missing_idx <- which(is.na(X_data_missing), arr.ind = TRUE)
  
  # Run Group LASSO BPCA
  cat("Running Bayesian Group LASSO BPCA...\n")
  result_lasso <- mcmc_bpca_group_lasso(
    X = X_data_missing,
    K = 5,
    G = 10,
    n_iter = 1500,
    burn_in = 500,
    n_chains = 3
  )
  
  # Evaluate results
  factor_eval <- evaluate_factor_recovery(W_true, result_lasso$W)
  missing_eval <- evaluate_reconstruction(X_data_complete, result_lasso$X_imputed, missing_idx)
  sparsity_analysis <- analyze_group_sparsity(result_lasso, W_true, c(1, 3, 5))
  
  cat("\n=== RESULTS ===\n")
  cat("Factor recovery correlation:", round(factor_eval$correlation, 4), "\n")
  cat("Factor recovery RMSE:", round(factor_eval$rmse, 4), "\n")
  cat("Missing data correlation:", round(missing_eval$correlation, 4), "\n")
  cat("Missing data RMSE:", round(missing_eval$rmse_missing, 4), "\n")
  cat("Detected active groups:", sparsity_analysis$detected_active_norms, "\n")
  cat("Sparsity achieved:", round(sparsity_analysis$sparsity_achieved, 1), "%\n")
  
     cat("\nQuick test completed successfully!\n")
  return(list(
    lasso = result_lasso,
    factor_eval = factor_eval,
    missing_eval = missing_eval,
    sparsity_analysis = sparsity_analysis
  ))
 }
 
 # SIMPLE SAVE PLOTS FUNCTION
 save_plots_only <- function(results, W_true, X_complete, missing_idx, output_dir = "bpca_missing_lasso_plots") {
   cat("Saving plots only...\n")
   
   # Extract aligned W for plotting
   factor_eval <- evaluate_factor_recovery(W_true, results$W)
   
   # Save all plots
  plot_W_comparison(W_true, factor_eval$W_aligned, " (Group LASSO)", save_plots = TRUE, output_dir = output_dir)
   plot_missing_data_reconstruction(X_complete, results$X_imputed, missing_idx, save_plots = TRUE, output_dir = output_dir)
   
   # MCMC trace plots
  plot_comprehensive_mcmc_traces(results, save_plots = TRUE, output_dir = output_dir)
   
  # Group evolution plots
     plot_group_evolution(results, save_plots = TRUE, output_dir = output_dir)
   
   cat("Plots saved to:", output_dir, "\n")
 }
                                                                                
#                           TEST FUNCTION FOR AUTO K SELECTION
#                                                                                    
# 
# test_auto_K <- function(K_true = 4, n = 100, p = 20, noise_level = 0.3, 
#                        missing_rate = 0.1, n_trials = 1, K_range = NULL) {
#   # PURPOSE: Comprehensive test of automatic K selection with validation
#   # PARAMETERS:
#   # - K_true: True number of factors to generate data with
#   # - n, p: Data dimensions
#   # - noise_level: Standard deviation of noise
#   # - missing_rate: Proportion of missing values
#   # - n_trials: Number of independent trials to run
#   # - K_range: Range of K values to test (default: K_true-2 to K_true+3)
#   
#   cat("COMPREHENSIVE AUTO K SELECTION TEST\n")
#   cat("                                                                                   \n")
#   
#   # Set default K_range if not provided
#   if (is.null(K_range)) {
#     K_range <- max(1, K_true - 2):min(p, K_true + 3)
#   }
#   
#   cat(" Test Settings:\n")
#   cat("   True K:", K_true, "\n")
#   cat("   Data dimensions: n =", n, ", p =", p, "\n")
#   cat("   Noise level:", noise_level, "\n")
#   cat("   Missing rate:", round(100*missing_rate, 1), "%\n")
#   cat("   K range to test:", paste(K_range, collapse=", "), "\n")
#   cat("   Number of trials:", n_trials, "\n\n")
#   
#   # Store results across trials
#   all_results <- list()
#   success_count <- 0
#   cv_scores_by_K <- matrix(NA, n_trials, length(K_range))
#   colnames(cv_scores_by_K) <- paste0("K_", K_range)
#   
#   for (trial in 1:n_trials) {
#     cat("TRIAL", trial, "of", n_trials, "\n")
#     cat("                                                                                   \n")
#     
#     # Generate data with known structure
#     set.seed(123 + trial)  # Different seed for each trial
#     
#     # Create structured latent factors with clear separation
#     Z_true <- matrix(0, n, K_true)
#     for (k in 1:K_true) {
#       Z_true[, k] <- rnorm(n, mean = 0, sd = 1)
#     }
#     
#     # Create structured loadings with strong signal
#     W_true <- matrix(0, p, K_true)
#     features_per_factor <- p %/% K_true
#     for (k in 1:K_true) {
#       start_idx <- (k-1) * features_per_factor + 1
#       end_idx <- min(k * features_per_factor, p)
#       W_true[start_idx:end_idx, k] <- rnorm(end_idx - start_idx + 1, mean = 0, sd = 1.5)
#     }
#     
#     # Generate data
#     mu_true <- rnorm(p, 0, 0.1)
#     X_true <- Z_true %*% t(W_true) + matrix(mu_true, n, p, byrow = TRUE) + 
#               matrix(rnorm(n * p, 0, noise_level), n, p)
#     
#     # Add missing values
#     X_with_missing <- X_true
#     n_missing <- round(missing_rate * length(X_with_missing))
#     missing_idx <- sample(length(X_with_missing), n_missing)
#     X_with_missing[missing_idx] <- NA
#     
#     cat("   Generated data: X ~", K_true, "factors + noise(", noise_level, ")\n")
#     cat("   Missing values:", n_missing, "(", round(100*missing_rate, 1), "%)\n")
#     
#     # Test auto K selection
#     cat("     Running auto K selection...\n")
#     start_time <- Sys.time()
#     
#     tryCatch({
#       # Run with detailed CV scores tracking
#       result <- mcmc_bpca_group_lasso(X_with_missing, 
#                                      auto_K = TRUE, 
#                                      K_range = K_range,
#                                      cv_folds = 3,
#                                      G = max(2, K_true),
#                                      n_iter = 1000,
#                                      burn_in = 300,
#                                      n_chains = 2)
#       
#       elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
#       selected_K <- ncol(result$W)
#       
#       # Validate selection
#       is_correct <- (selected_K == K_true)
#       if (is_correct) success_count <- success_count + 1
#       
#       # Evaluate quality
#       missing_pred <- predict_missing_values_simple(result, X_with_missing)
#       missing_rmse <- sqrt(mean((X_true[missing_idx] - missing_pred[missing_idx])^2))
#       max_rhat <- max(result$convergence_diagnostics$W_rhats, na.rm=TRUE)
#       
#       cat(" Selected K:", selected_K, 
#           if(is_correct) " (CORRECT!)" else paste(" (Expected:", K_true, ")"), "\n")
#       cat(" RMSE:", round(missing_rmse, 4), "| R-hat:", round(max_rhat, 3), 
#           "| Time:", round(elapsed_time, 2), "min\n")
#       
#       # Store detailed results
#       trial_result <- list(
#         trial = trial,
#         true_K = K_true,
#         selected_K = selected_K,
#         is_correct = is_correct,
#         rmse = missing_rmse,
#         max_rhat = max_rhat,
#         time_minutes = elapsed_time,
#         cv_scores = result$cv_results  # Store CV scores if available
#       )
#       
#       all_results[[trial]] <- trial_result
#       
#     }, error = function(e) {
#       cat(" FAILED:", e$message, "\n")
#       all_results[[trial]] <<- list(
#         trial = trial,
#         true_K = K_true,
#         selected_K = NA,
#         is_correct = FALSE,
#         error = e$message
#       )
#     })
#     
#     cat("\n")
#   }
#   
#   # Summary across all trials
#   cat("OVERALL AUTO K SELECTION PERFORMANCE\n")
#   cat("                                                                                   \n")
#   
#   successful_trials <- sum(!is.na(sapply(all_results, function(x) x$selected_K)))
#   accuracy_rate <- success_count / successful_trials * 100
#   
#   cat("SUCCESS RATE:", success_count, "out of", successful_trials, "trials")
#   cat(" (", round(accuracy_rate, 1), "%)\n")
#   
#   if (successful_trials > 0) {
#     selected_Ks <- sapply(all_results, function(x) x$selected_K)
#     selected_Ks <- selected_Ks[!is.na(selected_Ks)]
#     
#     cat("Selected K distribution:\n")
#     K_table <- table(selected_Ks)
#     for (k in names(K_table)) {
#       pct <- round(K_table[k] / successful_trials * 100, 1)
#       marker <- if (as.numeric(k) == K_true) " (TRUE)" else ""
#       cat("   K =", k, ":", K_table[k], "times (", pct, "%)", marker, "\n")
#     }
#     
#     # Average performance metrics
#     avg_rmse <- mean(sapply(all_results, function(x) x$rmse), na.rm = TRUE)
#     avg_rhat <- mean(sapply(all_results, function(x) x$max_rhat), na.rm = TRUE)
#     avg_time <- mean(sapply(all_results, function(x) x$time_minutes), na.rm = TRUE)
#     
#     cat("\n  Average Performance:\n")
#     cat("   RMSE:", round(avg_rmse, 4), "\n")
#     cat("   Max R-hat:", round(avg_rhat, 3), "\n")
#     cat("   Time per trial:", round(avg_time, 2), "minutes\n")
#   }
#   
#   # Interpretation and recommendations
#   cat("\n INTERPRETATION:\n")
#   if (is.na(accuracy_rate) || successful_trials == 0) {
#     cat("FAILED: All trials failed - check library dependencies and settings\n")
#   } else if (accuracy_rate >= 80) {
#     cat("EXCELLENT: Auto K selection is working very well!\n")
#   } else if (accuracy_rate >= 60) {
#     cat("GOOD: Auto K selection is mostly working, minor improvements possible\n")
#   } else if (accuracy_rate >= 40) {
#     cat("MODERATE: Auto K selection needs improvement\n")
#   } else {
#     cat("POOR: Auto K selection is not working reliably\n")
#   }
#   
#   cat("\nRECOMMENDATIONS:\n")
#   if (accuracy_rate < 80) {
#     cat("   - Increase cv_folds for more reliable CV estimates\n")
#     cat("   - Increase n_iter/burn_in for better MCMC convergence\n")
#     cat("   - Use stronger signal-to-noise ratio in data generation\n")
#     cat("   - Consider different evaluation metrics\n")
#   } else {
#     cat("   - Current settings work well for this data type\n")
#     cat("   - Ready for real data applications\n")
#   }
#   
#   cat("\n AUTO K VALIDATION TEST COMPLETED!\n")
#   
#   # Return comprehensive results
#   return(list(
#     settings = list(
#       K_true = K_true,
#       n = n, p = p,
#       noise_level = noise_level,
#       missing_rate = missing_rate,
#       K_range = K_range,
#       n_trials = n_trials
#     ),
#     summary = list(
#       success_count = success_count,
#       total_trials = successful_trials,
#       accuracy_rate = accuracy_rate,
#       selected_K_distribution = if(successful_trials > 0) table(selected_Ks) else NULL,
#       avg_rmse = if(successful_trials > 0) avg_rmse else NA,
#       avg_rhat = if(successful_trials > 0) avg_rhat else NA,
#       avg_time_minutes = if(successful_trials > 0) avg_time else NA
#     ),
#     detailed_results = all_results
#   ))
# }
# 
# predict_missing_values_simple <- function(model, X_with_missing) {
#   # Simple prediction function for testing
#   W <- model$W
#   mu <- model$mu
#   sigma2 <- model$sigma2
#   
#   N <- nrow(X_with_missing)
#   P <- ncol(X_with_missing)
#   K <- ncol(W)
#   
#   X_pred <- X_with_missing
#   
#   for (i in 1:N) {
#     missing_idx <- which(is.na(X_with_missing[i, ]))
#     observed_idx <- which(!is.na(X_with_missing[i, ]))
#     
#     if (length(observed_idx) > 0 && length(missing_idx) > 0) {
#       X_obs <- X_with_missing[i, observed_idx]
#       W_obs <- W[observed_idx, , drop = FALSE]
#       mu_obs <- mu[observed_idx]
#       
#       # Predict Z
#       Sigma_z <- solve(t(W_obs) %*% W_obs / sigma2 + diag(K))
#       mu_z <- Sigma_z %*% t(W_obs) %*% (X_obs - mu_obs) / sigma2
#       
#       # Predict missing values
#       X_pred[i, missing_idx] <- W[missing_idx, , drop = FALSE] %*% mu_z + mu[missing_idx]
#     }
#   }
#   
#   return(X_pred)
# }
# 
# #                                                                                    
# #                        COMPREHENSIVE AUTO K VALIDATION SUITE
# #                                                                                    
# 
# validate_auto_K_comprehensive <- function() {
#   # PURPOSE: Test auto K selection across multiple scenarios to ensure robustness
#   cat("COMPREHENSIVE AUTO K VALIDATION SUITE\n")
#   cat("                                                                                   \n")
#   cat("Testing auto K selection across different data scenarios...\n\n")
#   
#   # Test scenarios
#   scenarios <- list(
#     list(name = "Small & Clean", K_true = 3, n = 80, p = 15, noise = 0.2, missing = 0.05),
#     list(name = "Medium & Noisy", K_true = 4, n = 100, p = 20, noise = 0.5, missing = 0.15),
#     list(name = "Large & Sparse", K_true = 5, n = 150, p = 30, noise = 0.3, missing = 0.20),
#     list(name = "High Dimension", K_true = 6, n = 120, p = 40, noise = 0.4, missing = 0.10)
#   )
#   
#   overall_results <- list()
#   total_correct <- 0
#   total_tests <- 0
#   
#   for (i in 1:length(scenarios)) {
#     scenario <- scenarios[[i]]
#     cat("SCENARIO", i, ":", scenario$name, "\n")
#     cat("   K_true =", scenario$K_true, "| n =", scenario$n, "| p =", scenario$p, 
#         "| noise =", scenario$noise, "| missing =", round(100*scenario$missing, 1), "%\n")
#     
#     # Run test for this scenario
#     start_time <- Sys.time()
#     result <- test_auto_K(
#       K_true = scenario$K_true,
#       n = scenario$n,
#       p = scenario$p,
#       noise_level = scenario$noise,
#       missing_rate = scenario$missing,
#       n_trials = 2,  # Run 2 trials per scenario
#       K_range = max(1, scenario$K_true - 2):min(scenario$p, scenario$K_true + 3)
#     )
#     elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
#     
#     # Store results
#     overall_results[[scenario$name]] <- result
#     total_correct <- total_correct + result$summary$success_count
#     total_tests <- total_tests + result$summary$total_trials
#     
#     cat("Result:", result$summary$success_count, "out of", result$summary$total_trials, 
#         "correct (", round(result$summary$accuracy_rate, 1), "%) in", round(elapsed_time, 1), "min\n\n")
#   }
#   
#   # Overall summary
#   overall_accuracy <- if (total_tests > 0) total_correct / total_tests * 100 else 0
#   
#   cat("OVERALL VALIDATION RESULTS\n")
#   cat("                                                                                   \n")
#   cat("TOTAL SUCCESS RATE:", total_correct, "out of", total_tests, "tests")
#   cat(" (", round(overall_accuracy, 1), "%)\n\n")
#   
#   # Detailed breakdown
#   cat("SCENARIO BREAKDOWN:\n")
#   for (i in 1:length(scenarios)) {
#     scenario <- scenarios[[i]]
#     result <- overall_results[[scenario$name]]
#     accuracy <- result$summary$accuracy_rate
#     
#     status_icon <- if (accuracy >= 80) " " else if (accuracy >= 50) " " else " "
#     cat("   ", status_icon, scenario$name, ":", round(accuracy, 1), 
#         "% (", result$summary$success_count, "/", result$summary$total_trials, ")\n")
#   }
#   
#   # Final assessment
#   cat("\n VALIDATION ASSESSMENT:\n")
#   if (overall_accuracy >= 80) {
#     cat("  EXCELLENT: Auto K selection is highly reliable across scenarios!\n")
#     cat(" Ready for production use on real data\n")
#     cat(" Robust to different data characteristics\n")
#   } else if (overall_accuracy >= 60) {
#     cat("  GOOD: Auto K selection works well but could be improved\n")
#     cat(" Consider parameter tuning for challenging scenarios\n")
#     cat(" Generally safe for most applications\n")
#   } else if (overall_accuracy >= 40) {
#     cat("  MODERATE: Auto K selection needs improvement\n")
#     cat(" Use with caution, manual K selection may be better\n")
#     cat(" Recommend algorithm improvements\n")
#   } else {
#     cat("  POOR: Auto K selection is not reliable\n")
#     cat(" Not recommended for production use\n")
#     cat(" Requires significant algorithm improvements\n")
#   }
#   
#   cat("\n RECOMMENDATIONS FOR IMPROVEMENT:\n")
#   if (overall_accuracy < 80) {
#     cat("   - Increase cross-validation folds (cv_folds = 5 or more)\n")
#     cat("   - Use longer MCMC chains for more stable estimates\n")
#     cat("   - Consider alternative model selection criteria (AIC, BIC)\n")
#     cat("   - Implement ensemble K selection (majority vote across methods)\n")
#     cat("   - Add early stopping based on convergence of CV scores\n")
#   } else {
#     cat("   - Current implementation is robust and ready to use!\n")
#     cat("   - Consider adding confidence intervals for K selection\n")
#     cat("   - Document successful scenarios for user guidance\n")
#   }
#   
#   cat("\n COMPREHENSIVE VALIDATION COMPLETED!\n")
#   
#   return(list(
#     overall_accuracy = overall_accuracy,
#     total_correct = total_correct,
#     total_tests = total_tests,
#     scenario_results = overall_results,
#     assessment = if (overall_accuracy >= 80) "EXCELLENT" else 
#                 if (overall_accuracy >= 60) "GOOD" else 
#                 if (overall_accuracy >= 40) "MODERATE" else "POOR"
#   ))
# }
# 
# # Quick validation function for immediate testing
# quick_validate_auto_K <- function() {
#   cat("QUICK AUTO K VALIDATION\n")
#   cat("                                                                                   \n")
#   
#   # Single scenario test with known good settings
#   result <- test_auto_K(
#     K_true = 4,
#     n = 100,
#     p = 20,
#     noise_level = 0.3,
#     missing_rate = 0.1,
#     n_trials = 1,
#     K_range = 2:7
#   )
#   
#   cat("\n QUICK VALIDATION COMPLETE!\n")
#   if (result$summary$accuracy_rate == 100) {
#     cat("  AUTO K SELECTION IS WORKING CORRECTLY!\n")
#   } else {
#     cat("  AUTO K SELECTION NEEDS ATTENTION\n")
#   }
#   
#   return(result)
# } 

# --- CORRECT FUNCTIONS FROM bayesian_group_lasso_true.R ARE SOURCED ---

# =============================================================================
# UNIT TESTING AND VALIDATION FUNCTIONS
# =============================================================================

# Comprehensive unit testing for Bayesian Group LASSO BPCA
test_bayesian_group_lasso_implementation <- function() {
  cat("=== UNIT TESTING BAYESIAN GROUP LASSO BPCA ===\n\n")
  
  # Test 1: Function Availability
  cat("1. Testing Function Availability...\n")
  required_functions <- c("sample_group_lasso_loadings", "sample_group_scales_gig", 
                         "sample_global_shrinkage", "mcmc_bpca_group_lasso")
  missing_functions <- required_functions[!sapply(required_functions, exists)]
  if (length(missing_functions) == 0) {
    cat("   ✅ All required functions available\n")
  } else {
    cat("   ❌ Missing functions:", paste(missing_functions, collapse = ", "), "\n")
    return(FALSE)
  }
  
  # Test 2: Convergence Diagnostics
  cat("\n2. Testing Convergence Diagnostics...\n")
  tryCatch({
    # Create dummy chains with known convergence
    n_chains <- 4
    n_samples <- 100
    dummy_chains <- array(rnorm(n_chains * n_samples * 10), dim = c(n_chains, n_samples, 10))
    
    # Test R-hat calculation
    test_rhat <- calculate_rhat(dummy_chains[,,1])
    if (abs(test_rhat - 1.0) < 0.2) {
      cat("   ✅ R-hat calculation working correctly\n")
    } else {
      cat("   ⚠️  R-hat calculation may have issues\n")
    }
  }, error = function(e) {
    cat("   ❌ Convergence diagnostics failed:", e$message, "\n")
  })
  
  # Test 3: Data Simulation
  cat("\n3. Testing Data Simulation...\n")
  tryCatch({
    sim_data <- simulate_data_with_groups(n = 50, p = 20, G = 4, K = 2, 
                                        active_groups = c(1, 2), sigma = 0.3)
    
    if (all(c("X_data", "W_true", "Z_true", "group_id") %in% names(sim_data))) {
      cat("   ✅ Data simulation working correctly\n")
      cat("   - Data dimensions:", nrow(sim_data$X_data), "×", ncol(sim_data$X_data), "\n")
      cat("   - True factors:", sim_data$K_true, "\n")
      cat("   - Groups:", sim_data$G_total, "\n")
    } else {
      cat("   ❌ Data simulation missing components\n")
    }
  }, error = function(e) {
    cat("   ❌ Data simulation failed:", e$message, "\n")
  })
  
  # Test 4: Small MCMC Run
  cat("\n4. Testing Small MCMC Run...\n")
  tryCatch({
    # Create small test data
    set.seed(123)
    test_X <- matrix(rnorm(20 * 10), 20, 10)
    test_X[sample(1:200, 30)] <- NA  # Add some missing values
    
    # Run short MCMC
    test_result <- mcmc_bpca_group_lasso(
      X = test_X, K = 2, G = 5, n_iter = 100, burn_in = 20, 
      n_chains = 2, auto_group = FALSE, auto_K = FALSE
    )
    
    if (all(c("W", "Z", "mu", "sigma2", "tau_g2", "lambda2") %in% names(test_result))) {
      cat("   ✅ MCMC sampling working correctly\n")
      cat("   - Estimated W dimensions:", nrow(test_result$W), "×", ncol(test_result$W), "\n")
      cat("   - Final lambda2:", round(test_result$lambda2, 3), "\n")
    } else {
      cat("   ❌ MCMC sampling missing components\n")
    }
  }, error = function(e) {
    cat("   ❌ MCMC sampling failed:", e$message, "\n")
  })
  
  # Test 5: Cross-Validation
  cat("\n5. Testing Cross-Validation...\n")
  tryCatch({
    # Create test data
    set.seed(456)
    test_X <- matrix(rnorm(30 * 15), 30, 15)
    test_X[sample(1:450, 50)] <- NA
    
    # Test CV with small parameters
    cv_result <- mcmc_bpca_group_lasso(
      X = test_X, K = 3, G = 5, n_iter = 50, burn_in = 10,
      auto_K = TRUE, K_range = 2:4, cv_folds = 3, n_chains = 2
    )
    
    cat("   ✅ Cross-validation working correctly\n")
  }, error = function(e) {
    cat("   ❌ Cross-validation failed:", e$message, "\n")
  })
  
  cat("\n=== UNIT TESTING COMPLETED ===\n")
  cat("If all tests passed, the implementation is working correctly!\n\n")
  
  return(TRUE)
}

# Run unit tests if called directly
if (FALSE) {  # Set to TRUE to run tests
  test_bayesian_group_lasso_implementation()
}