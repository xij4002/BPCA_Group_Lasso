library(mvtnorm)
library(GIGrvg)
source("bayesian_group_lasso_true copy.R")

# Setup test parameters
N <- 50
P <- 10
K <- 3
G <- 2
group_id <- rep(1:G, each = P / G)

set.seed(42)
Z <- matrix(rnorm(N * K), nrow = N, ncol = K)
X_residual <- matrix(rnorm(N * P), nrow = N, ncol = P)
sigma2 <- 1.5
tau_g_squared <- c(0.5, 1.2)

# Test group 1
g <- 1
vars_in_group_g1 <- which(group_id == g)
d_g_g1 <- length(vars_in_group_g1)
X_residual_g1 <- X_residual[, vars_in_group_g1, drop = FALSE]

# Calculate theoretical values
Z_T_Z <- crossprod(Z)
Omega_g1_theory <- Z_T_Z + (1 / tau_g_squared[g]) * diag(K)
B_g1 <- crossprod(X_residual_g1, Z)
Omega_g1_inv_theory <- solve(Omega_g1_theory)
mu_W_g1_theory <- B_g1 %*% Omega_g1_inv_theory
Sigma_W_g1_theory <- sigma2 * kronecker(Omega_g1_inv_theory, diag(d_g_g1))

mu_1_1_theory <- mu_W_g1_theory[1, 1]
var_1_1_theory <- Sigma_W_g1_theory[1, 1]
sd_1_1_theory <- sqrt(var_1_1_theory)

# Monte Carlo sampling for W[1,1]
n_samples <- 20000
w11_samples <- numeric(n_samples)

for (i in 1:n_samples) {
  W_sample <- sample_group_lasso_loadings(X_residual, Z, sigma2, tau_g_squared, group_id, G)
  w11_samples[i] <- W_sample[1, 1]
}

# Compare results
cat("Theoretical mean (W[1,1]):", mu_1_1_theory, "\n")
cat("Empirical mean (W[1,1]):", mean(w11_samples), "\n\n")
cat("Theoretical std dev (W[1,1]):", sd_1_1_theory, "\n")
cat("Empirical std dev (W[1,1]):", sd(w11_samples), "\n\n")

# Plot W[1,1] distribution
hist(w11_samples, breaks = 50, freq = FALSE,
     main = "Empirical vs Theoretical Distribution of W[1,1]",
     xlab = "W[1,1] samples")
curve(dnorm(x, mean = mu_1_1_theory, sd = sd_1_1_theory),
      from = min(w11_samples), to = max(w11_samples),
      col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Empirical (histogram)", "Theoretical (red line)"),
       col = c("black", "red"), lty = 1, lwd = c(1, 2))

cat("--- W[1,1] test plot generated ---\n")
cat("--- Press Enter to continue testing sample_group_scales_gig ---\n")
invisible(readline())

# Test sample_group_scales_gig
cat("\n--- Starting sample_group_scales_gig test ---\n")

set.seed(43)
W_fixed <- matrix(rnorm(P * K), nrow = P, ncol = K)
sigma2_fixed <- 1.5
lambda2_fixed <- 2.0

g <- 1
vars_in_group_g1 <- which(group_id == g)
W_g1_fixed <- W_fixed[vars_in_group_g1, , drop = FALSE]

# Theoretical GIG parameters
p_theory <- 1/2
a_theory <- lambda2_fixed
b_theory <- sum(W_g1_fixed^2) / sigma2_fixed

# Monte Carlo sampling for tau_g^2[1]
n_samples_gig <- 20000
tau_g1_samples <- numeric(n_samples_gig)

for (i in 1:n_samples_gig) {
  tau_g_all <- sample_group_scales_gig(W_fixed, group_id, sigma2_fixed, lambda2_fixed, G, K)
  tau_g1_samples[i] <- tau_g_all[1]
}

# Plot tau_g^2[1] distribution
hist(tau_g1_samples, breaks = 50, freq = FALSE,
     main = "Empirical vs Theoretical Distribution of tau_g^2[1]",
     xlab = "tau_g^2[1] samples",
     xlim = c(0, quantile(tau_g1_samples, 0.99)))
curve(GIGrvg::dgig(x, lambda = p_theory, chi = b_theory, psi = a_theory),
      from = 0.001, to = max(tau_g1_samples),
      col = "blue", lwd = 2, add = TRUE, n = 1001)
legend("topright", legend = c("Empirical (histogram)", "Theoretical (blue line)"),
       col = c("black", "blue"), lty = 1, lwd = c(1, 2))

cat("--- tau_g^2[1] test plot generated ---\n")
cat("--- All tests completed ---\n")