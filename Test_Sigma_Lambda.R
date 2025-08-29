library(invgamma)
library(ggplot2)

set.seed(123)

sample_global_shrinkage <- function(tau_g_squared, group_sizes, lambda_a, lambda_b, K) {
  G <- length(tau_g_squared)
  
  shape_increment <- 0
  for (g in 1:G) {
    m_g <- group_sizes[g] * K
    shape_increment <- shape_increment + (m_g + 1) / 2
  }
  
  shape_post <- lambda_a + shape_increment
  rate_post <- lambda_b + 0.5 * sum(tau_g_squared)
  
  lambda2 <- rgamma(1, shape = shape_post, rate = rate_post)
  
  if (!is.finite(lambda2) || lambda2 <= 1e-9) {
    lambda2 <- 0.1
  }
  
  return(lambda2)
}

sample_residual_variance <- function(X_residual, W, Z, sigma_a, sigma_b, tau_g_squared, group_id) {
  N <- nrow(X_residual)
  P <- ncol(X_residual)
  K <- ncol(W)
  
  X_fitted <- Z %*% t(W)
  residual_matrix <- X_residual - X_fitted
  data_sse <- sum(residual_matrix^2)
  
  prior_penalty <- 0
  total_params <- 0
  G <- length(tau_g_squared)
  for (g in 1:G) {
    vars_in_group <- which(group_id == g)
    if (length(vars_in_group) > 0) {
      W_g <- W[vars_in_group, , drop = FALSE]
      prior_penalty <- prior_penalty + sum(W_g^2) / tau_g_squared[g]
      total_params <- total_params + length(vars_in_group) * K
    }
  }
  
  shape_post <- sigma_a + (N * P + total_params) / 2
  scale_post <- sigma_b + (data_sse + prior_penalty) / 2
  sigma2 <- rinvgamma(1, shape = shape_post, scale = scale_post)
  # sigma2 <- 1 / rgamma(1, shape = shape_post, rate = scale_post)
  
  if (!is.finite(sigma2) || sigma2 <= 1e-9) {
    sigma2 <- 0.01
  }
  
  return(sigma2)
}

cat("=== Parameter Setup ===\n")
N <- 200  # Number of samples
P <- 30   # Number of variables  
K <- 3    # Number of factors
G <- 30   # Number of groups (one variable per group)

# Prior parameters
sigma_a_prior <- 1.0
sigma_b_prior <- 1.0
lambda_a_prior <- 1.0
lambda_b_prior <- 1.0

# True parameters
sigma2_true <- 2.0
lambda2_true <- 1.5
tau_g_squared_true <- rgamma(G, shape = 2, rate = lambda2_true)  # Shrinkage parameters for each group
group_id <- 1:G  # Each variable belongs to one group

cat("N =", N, ", P =", P, ", K =", K, ", G =", G, "\n")
cat("sigma2_true =", sigma2_true, "\n")
cat("lambda2_true =", lambda2_true, "\n")
cat("tau_g_squared_true (first 5):", head(tau_g_squared_true), "...\n\n")

# ===================================================================
# Generate true data
# ===================================================================
cat("=== Generating Simulated Data ===\n")

# Generate W matrix (loading matrix)
W_true <- matrix(0, P, K)
for (g in 1:G) {
  vars_in_group <- which(group_id == g)
  for (k in 1:K) {
    W_true[vars_in_group, k] <- rnorm(length(vars_in_group), 
                                      mean = 0, 
                                      sd = sqrt(sigma2_true * tau_g_squared_true[g]))
  }
}

# Generate Z matrix (factor scores)
Z_true <- matrix(rnorm(N * K, 0, 1), N, K)

# Generate observed data
X_fitted_true <- Z_true %*% t(W_true)
noise <- matrix(rnorm(N * P, 0, sqrt(sigma2_true)), N, P)
X_true <- X_fitted_true + noise

cat("Data generation completed\n")
cat("X_true dimensions:", dim(X_true), "\n")
cat("W_true dimensions:", dim(W_true), "\n")
cat("Z_true dimensions:", dim(Z_true), "\n\n")

# ===================================================================
# Test 1: sample_residual_variance
# ===================================================================
cat("=== Test 1: sample_residual_variance ===\n")

# Calculate theoretical posterior parameters
X_fitted_theory <- Z_true %*% t(W_true)
residual_matrix_theory <- X_true - X_fitted_theory
data_sse_theory <- sum(residual_matrix_theory^2)

prior_penalty_theory <- 0
total_params_theory <- 0
for (g in 1:G) {
  vars_in_group <- which(group_id == g)
  W_g <- W_true[vars_in_group, , drop = FALSE]
  prior_penalty_theory <- prior_penalty_theory + sum(W_g^2) / tau_g_squared_true[g]
  total_params_theory <- total_params_theory + length(vars_in_group) * K
}

shape_theory <- sigma_a_prior + (N * P + total_params_theory) / 2
scale_theory <- sigma_b_prior + (data_sse_theory + prior_penalty_theory) / 2

cat("Theoretical posterior parameters:\n")
cat("  shape =", shape_theory, "\n")
cat("  scale =", scale_theory, "\n")
cat("  theoretical mean =", scale_theory / (shape_theory - 1), "\n\n")

# Generate samples from function
n_samples <- 10000
sigma2_samples <- numeric(n_samples)

cat("Generating", n_samples, "samples...\n")
for (i in 1:n_samples) {
  sigma2_samples[i] <- sample_residual_variance(
    X_residual = X_true,
    W = W_true,
    Z = Z_true,
    sigma_a = sigma_a_prior,
    sigma_b = sigma_b_prior,
    tau_g_squared = tau_g_squared_true,
    group_id = group_id
  )
}

# Generate samples from theoretical distribution
sigma2_theory_samples <- rinvgamma(n_samples, shape = shape_theory, scale = scale_theory)

# Statistical comparison
cat("Sample statistics:\n")
cat("  function sample mean =", mean(sigma2_samples), "\n")
cat("  theoretical sample mean =", mean(sigma2_theory_samples), "\n")
cat("  function sample variance =", var(sigma2_samples), "\n")
cat("  theoretical sample variance =", var(sigma2_theory_samples), "\n\n")

# KS test
ks_test1 <- ks.test(sigma2_samples, sigma2_theory_samples)
cat("KS test results:\n")
cat("  D statistic =", ks_test1$statistic, "\n")
cat("  p-value =", ks_test1$p.value, "\n")
cat("  conclusion:", ifelse(ks_test1$p.value > 0.05, "✓ Function correct (p > 0.05)", "✗ Function may have issues (p ≤ 0.05)"), "\n\n")

# ===================================================================
# Test 2: sample_global_shrinkage  
# ===================================================================
cat("=== Test 2: sample_global_shrinkage ===\n")

# Calculate theoretical posterior parameters
group_sizes <- rep(1, G)  # One variable per group
shape_increment_theory <- 0
for (g in 1:G) {
  m_g <- group_sizes[g] * K  # Number of parameters per group
  shape_increment_theory <- shape_increment_theory + (m_g + 1) / 2
}

shape_theory2 <- lambda_a_prior + shape_increment_theory
rate_theory2 <- lambda_b_prior + 0.5 * sum(tau_g_squared_true)

cat("Theoretical posterior parameters:\n")
cat("  shape =", shape_theory2, "\n")
cat("  rate =", rate_theory2, "\n")
cat("  theoretical mean =", shape_theory2 / rate_theory2, "\n\n")

# Generate samples from function
lambda2_samples <- numeric(n_samples)

cat("Generating", n_samples, "samples...\n")
for (i in 1:n_samples) {
  lambda2_samples[i] <- sample_global_shrinkage(
    tau_g_squared = tau_g_squared_true,
    group_sizes = group_sizes,
    lambda_a = lambda_a_prior,
    lambda_b = lambda_b_prior,
    K = K
  )
}

# Generate samples from theoretical distribution
lambda2_theory_samples <- rgamma(n_samples, shape = shape_theory2, rate = rate_theory2)

# Statistical comparison
cat("Sample statistics:\n")
cat("  function sample mean =", mean(lambda2_samples), "\n")
cat("  theoretical sample mean =", mean(lambda2_theory_samples), "\n")
cat("  function sample variance =", var(lambda2_samples), "\n")
cat("  theoretical sample variance =", var(lambda2_theory_samples), "\n\n")

# KS test
ks_test2 <- ks.test(lambda2_samples, lambda2_theory_samples)
cat("KS test results:\n")
cat("  D statistic =", ks_test2$statistic, "\n")
cat("  p-value =", ks_test2$p.value, "\n")
cat("  conclusion:", ifelse(ks_test2$p.value > 0.05, "✓ Function correct (p > 0.05)", "✗ Function may have issues (p ≤ 0.05)"), "\n\n")

# ===================================================================
# Visualization comparison
# ===================================================================
cat("=== Generating Visualization ===\n")

# Create comparison data frames
comparison_df1 <- data.frame(
  samples = c(sigma2_samples, sigma2_theory_samples),
  source = rep(c("Function Output", "Theoretical Distribution"), each = n_samples),
  test = "sample_residual_variance"
)

comparison_df2 <- data.frame(
  samples = c(lambda2_samples, lambda2_theory_samples),
  source = rep(c("Function Output", "Theoretical Distribution"), each = n_samples),
  test = "sample_global_shrinkage"
)

# Plot sigma2 comparison
g1 <- ggplot(comparison_df1, aes(x = samples, fill = source)) +
  geom_density(alpha = 0.6) +
  labs(title = "sample_residual_variance Validation",
       subtitle = paste("KS p-value =", round(ks_test1$p.value, 4)),
       x = "sigma2",
       y = "Density",
       fill = "Sample Source") +
  theme_minimal()

# Plot lambda2 comparison
g2 <- ggplot(comparison_df2, aes(x = samples, fill = source)) +
  geom_density(alpha = 0.6) +
  labs(title = "sample_global_shrinkage Validation",
       subtitle = paste("KS p-value =", round(ks_test2$p.value, 4)),
       x = "lambda2",
       y = "Density",
       fill = "Sample Source") +
  theme_minimal()

print(g1)
print(g2)

# ===================================================================
# Final summary
# ===================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("FINAL TEST RESULTS:\n")
cat("1. sample_residual_variance:", 
    ifelse(ks_test1$p.value > 0.05, "✓ PASSED", "✗ FAILED"), 
    "(p =", round(ks_test1$p.value, 4), ")\n")
cat("2. sample_global_shrinkage:", 
    ifelse(ks_test2$p.value > 0.05, "✓ PASSED", "✗ FAILED"), 
    "(p =", round(ks_test2$p.value, 4), ")\n")

if (ks_test1$p.value > 0.05 && ks_test2$p.value > 0.05) {
  cat("\nCONGRATULATIONS! Both functions correctly implement their respective posterior sampling.\n")
} else {
  cat("\n⚠️  Functions with issues need further investigation of implementation details.\n")
}

cat("Test parameters: N =", N, ", P =", P, ", K =", K, ", G =", G, "\n")
cat("Prior parameters: sigma_a = sigma_b = lambda_a = lambda_b = 1.0\n")