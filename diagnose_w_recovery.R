# =============================================================================
# DIAGNOSE W RECOVERY ISSUES
# =============================================================================

# Load required libraries and functions
source("bpca_missing_lasso.R")

# Set seed for reproducibility
set.seed(123)

# =============================================================================
# STEP 1: GENERATE DATA AND ANALYZE TRUE W STRUCTURE
# =============================================================================

cat("=== DIAGNOSING W RECOVERY ISSUES ===\n\n")

# Generate data with K=3
cat("1. Analyzing K=3 data structure...\n")
sim_data_k3 <- simulate_data_with_groups(
  n = 200,      # 200 subjects
  p = 100,      # 100 variables
  G = 100,      # 100 groups
  K = 3,        # 3 factors
  active_groups = c(1, 2, 3),
  sigma = 0.3
)

W_true_k3 <- sim_data_k3$W_true
cat("   True W (K=3) dimensions:", dim(W_true_k3), "\n")
cat("   True W (K=3) non-zero elements:", sum(W_true_k3 != 0), "/", length(W_true_k3), "\n")
cat("   True W (K=3) sparsity:", round(100 * (1 - sum(W_true_k3 != 0) / length(W_true_k3)), 1), "%\n")
cat("   True W (K=3) Frobenius norm:", round(norm(W_true_k3, "F"), 4), "\n")
cat("   True W (K=3) column norms:", round(apply(W_true_k3, 2, function(x) sqrt(sum(x^2))), 4), "\n")
cat("   True W (K=3) range:", round(range(W_true_k3), 4), "\n\n")

# Generate data with K=5
cat("2. Analyzing K=5 data structure...\n")
sim_data_k5 <- simulate_data_with_groups(
  n = 200,      # 200 subjects
  p = 100,      # 100 variables
  G = 100,      # 100 groups
  K = 5,        # 5 factors
  active_groups = c(1, 2, 3, 4, 5),
  sigma = 0.3
)

W_true_k5 <- sim_data_k5$W_true
cat("   True W (K=5) dimensions:", dim(W_true_k5), "\n")
cat("   True W (K=5) non-zero elements:", sum(W_true_k5 != 0), "/", length(W_true_k5), "\n")
cat("   True W (K=5) sparsity:", round(100 * (1 - sum(W_true_k5 != 0) / length(W_true_k5)), 1), "%\n")
cat("   True W (K=5) Frobenius norm:", round(norm(W_true_k5, "F"), 4), "\n")
cat("   True W (K=5) column norms:", round(apply(W_true_k5, 2, function(x) sqrt(sum(x^2))), 4), "\n")
cat("   True W (K=5) range:", round(range(W_true_k5), 4), "\n\n")

# =============================================================================
# STEP 2: TEST DIFFERENT PRIOR SPECIFICATIONS
# =============================================================================

cat("3. Testing different prior specifications for K=3...\n")

# Test different lambda_a and lambda_b combinations
prior_configs <- list(
  list(name = "Strong LASSO", lambda_a = 2.0, lambda_b = 1.0),
  list(name = "Medium LASSO", lambda_a = 1.0, lambda_b = 0.5),
  list(name = "Weak LASSO", lambda_a = 0.5, lambda_b = 0.25),
  list(name = "Very Weak LASSO", lambda_a = 0.1, lambda_b = 0.05)
)

# Introduce missing values for K=3
X_data_k3 <- sim_data_k3$X_data
set.seed(456)
missing_rate <- 0.10
n_missing <- round(nrow(X_data_k3) * ncol(X_data_k3) * missing_rate)
missing_indices <- sample(1:length(X_data_k3), n_missing)
X_data_k3_missing <- X_data_k3
X_data_k3_missing[missing_indices] <- NA

# Test each prior configuration
for (config in prior_configs) {
  cat("   Testing", config$name, "...\n")
  
  # Run MCMC with shorter iterations for testing
  result <- mcmc_bpca_group_lasso(
    X = X_data_k3_missing,
    K = 3,
    G = 100,
    n_iter = 1000,  # Shorter for testing
    burn_in = 250,
    n_chains = 2,   # Fewer chains for speed
    lambda_a = config$lambda_a,
    lambda_b = config$lambda_b,
    auto_group = FALSE,
    auto_K = FALSE
  )
  
  # Evaluate W recovery
  factor_eval <- evaluate_factor_recovery(W_true_k3, result$W)
  
  cat("     - Factor recovery correlation:", round(factor_eval$correlation, 4), "\n")
  cat("     - Factor recovery RMSE:", round(factor_eval$rmse, 4), "\n")
  cat("     - Final lambda2:", round(result$lambda2, 3), "\n")
  cat("     - Estimated W Frobenius norm:", round(norm(result$W, "F"), 4), "\n")
  cat("     - Norm ratio (est/true):", round(norm(result$W, "F") / norm(W_true_k3, "F"), 4), "\n\n")
}

cat("4. Testing different prior specifications for K=5...\n")

# Introduce missing values for K=5
X_data_k5 <- sim_data_k5$X_data
set.seed(456)
n_missing <- round(nrow(X_data_k5) * ncol(X_data_k5) * missing_rate)
missing_indices <- sample(1:length(X_data_k5), n_missing)
X_data_k5_missing <- X_data_k5
X_data_k5_missing[missing_indices] <- NA

# Test each prior configuration for K=5
for (config in prior_configs) {
  cat("   Testing", config$name, "...\n")
  
  # Run MCMC with shorter iterations for testing
  result <- mcmc_bpca_group_lasso(
    X = X_data_k5_missing,
    K = 5,
    G = 100,
    n_iter = 1000,  # Shorter for testing
    burn_in = 250,
    n_chains = 2,   # Fewer chains for speed
    lambda_a = config$lambda_a,
    lambda_b = config$lambda_b,
    auto_group = FALSE,
    auto_K = FALSE
  )
  
  # Evaluate W recovery
  factor_eval <- evaluate_factor_recovery(W_true_k5, result$W)
  
  cat("     - Factor recovery correlation:", round(factor_eval$correlation, 4), "\n")
  cat("     - Factor recovery RMSE:", round(factor_eval$rmse, 4), "\n")
  cat("     - Final lambda2:", round(result$lambda2, 3), "\n")
  cat("     - Estimated W Frobenius norm:", round(norm(result$W, "F"), 4), "\n")
  cat("     - Norm ratio (est/true):", round(norm(result$W, "F") / norm(W_true_k5, "F"), 4), "\n\n")
}

cat("=== DIAGNOSIS COMPLETED ===\n")
cat("Check the results above to find the best prior configuration.\n") 