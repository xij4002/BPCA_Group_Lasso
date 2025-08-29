# source("bayesian_group_lasso_true.R")
source("bayesian_group_lasso_true copy.R")
library(mvtnorm)
library(tidyverse)
library(pheatmap)
library(gridExtra)
library(grid)
library(vegan)
library(coda)

library(Matrix)
library(mvtnorm)

# simulate_data_with_groups: Takes dimensions and a sparsity map as input, generates a complete simulated dataset with true W, Z, and mu.
# sample_z: Takes W, X, mu, and sigma2 as input, performs an efficient, vectorized Gibbs sample for the entire latent factor matrix Z, and outputs the new Z.
# update_W_Z_blocks: Takes W, Z, and other parameters as input, performs a block Gibbs sample for a single, randomly chosen factor pair (W[,k], Z[,k]), and outputs the updated W and Z.
# mcmc_bpca_group_lasso: The main MCMC engine. It takes data and model parameters as input, runs multiple MCMC chains using a mix of vectorized and block Gibbs sampling, and outputs the raw posterior samples.
# align_factor_chains: Takes the raw MCMC output for W and Z from multiple chains as input, aligns them using Procrustes analysis on the posterior means, and outputs the aligned sample chains.
# calculate_convergence_diagnostics: Takes aligned MCMC chains as input, calculates R-hat convergence statistics for key parameters, and outputs a summary.
# evaluate_factor_recovery: Takes the true W and the estimated posterior mean W as input, aligns them, calculates correlation and RMSE, and outputs recovery metrics.
# plot_W_comparison: Takes true and estimated W matrices as input, performs normalization tests, generates comparison heatmaps, and outputs the plots.
# debug_scale_correction: A post-hoc utility that takes posterior chains as input, corrects for scale ambiguity by normalizing Z's variance, and outputs the corrected W.
# diagnosis_report: A top-level wrapper that takes MCMC results and true parameters as input, runs all diagnostic, recovery, and plotting functions, and outputs a comprehensive summary report.

simulate_data_with_groups <- function(n, p, G, K, active_map, sigma, mu_sd = 1, w_min = 0.6, w_max = 1.0) {
  if (p %% G != 0) stop("p must be divisible by G")
  active_factors <- as.integer(names(active_map))
  if (any(active_factors > K)) stop("Active factor index cannot exceed K")
  if (any(unlist(active_map) > G)) stop("Active group index cannot exceed G")
  
  group_size <- p / G
  group_id <- rep(1:G, each = group_size)
  
  W_true <- matrix(0, nrow = p, ncol = K)
  
  for (k in active_factors) {
    associated_groups <- active_map[[as.character(k)]]
    for (g in associated_groups) {
      rows_in_group <- which(group_id == g)
      sign <- sample(c(-1, 1), 1)
      W_true[rows_in_group, k] <- runif(length(rows_in_group), w_min, w_max) * sign
    }
  }
  
  Z_true <- matrix(rnorm(n * K), nrow = n, ncol = K)
  mu_true <- rnorm(p, 0, mu_sd)
  noise <- matrix(rnorm(n * p, 0, sigma), nrow = n, ncol = p)
  
  X_data_complete <- Z_true %*% t(W_true) + matrix(rep(mu_true, n), nrow = n, byrow = TRUE) + noise
  
  return(list(
    X_data = X_data_complete,
    W_true = W_true,
    Z_true = Z_true,
    group_id = group_id,
    K_true = K,
    G_total = G
  ))
}

library(Matrix)
library(mvtnorm)

sample_z <- function(W, X_work, mu, sigma2, N, K) {
  WtW <- t(W) %*% W
  mu_matrix <- matrix(rep(mu, N), nrow = N, byrow = TRUE)
  X_centered <- X_work - mu_matrix
  
  precision_Z <- WtW / sigma2 + diag(K)
  
  min_eig <- min(eigen(precision_Z, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig < 1e-10) {
    precision_Z <- precision_Z + (1e-8 - min_eig) * diag(K)
  }
  
  Sigma_Z <- solve(precision_Z)
  mu_Z <- (X_centered %*% W) %*% Sigma_Z / sigma2
  
  L <- chol(Sigma_Z)
  noise <- matrix(rnorm(N * K), N, K)
  Z_samples <- mu_Z + noise %*% L
  
  return(Z_samples)
}

update_W_Z_blocks <- function(X_work, W, Z, mu, sigma2, group_id, G, K, no_shrinkage) {
  N <- nrow(X_work)
  P <- ncol(X_work)
  k <- sample(1:K, 1)
  
  if (K == 1) {
    X_residual <- X_work - matrix(rep(mu, each = N), nrow = N)
  } else {
    other_factors <- setdiff(1:K, k)
    X_residual <- X_work - matrix(rep(mu, each = N), nrow = N) - 
      Z[, other_factors, drop = FALSE] %*% t(W[, other_factors, drop = FALSE])
  }
  
  ZkZk <- sum(Z[, k]^2)
  if (no_shrinkage) {
    precision_w <- ZkZk / sigma2 + 1e-6
  } else {
    precision_w <- ZkZk / sigma2 + 1
  }
  Sigma_w <- 1 / precision_w
  
  mu_w <- (t(X_residual) %*% Z[, k]) * Sigma_w / sigma2
  W[, k] <- mu_w + rnorm(P, 0, sqrt(Sigma_w))
  
  if (K == 1) {
    X_residual <- X_work - matrix(rep(mu, each = N), nrow = N)
  } else {
    other_factors <- setdiff(1:K, k)
    X_residual <- X_work - matrix(rep(mu, each = N), nrow = N) - 
      Z[, other_factors, drop = FALSE] %*% t(W[, other_factors, drop = FALSE])
  }
  
  WkWk <- sum(W[, k]^2)
  precision_z <- WkWk / sigma2 + 1
  Sigma_z <- 1 / precision_z
  
  mu_z <- (X_residual %*% W[, k]) * Sigma_z / sigma2
  Z[, k] <- mu_z + rnorm(N, 0, sqrt(Sigma_z))
  
  return(list(W = W, Z = Z, updated_factor = k))
}

mcmc_bpca_group_lasso <- function(X, K, G, n_iter, burn_in, n_chains = 3, 
                                  mu_prior_var = 1000, 
                                  sigma_a = 0.1, sigma_b = 0.1,
                                  lambda_a = 0.1, lambda_b = 0.1,
                                  fix_sigma2 = NA,
                                  no_shrinkage = FALSE,
                                  use_block_sampling = TRUE) {
  N <- nrow(X)
  P <- ncol(X)
  
  if (P %% G != 0) {
    stop("P must be divisible by G for equal group sizes")
  }
  
  group_size <- P %/% G
  group_id <- rep(1:G, each = group_size)
  group_sizes_vec <- rep(group_size, G)
  
  missing_idx <- which(is.na(X), arr.ind = TRUE)
  X_work <- X
  if (nrow(missing_idx) > 0) {
    for (j in 1:P) {
      missing_j <- is.na(X[, j])
      if (any(missing_j)) {
        col_mean <- mean(X[, j], na.rm = TRUE)
        X_work[missing_j, j] <- ifelse(is.finite(col_mean), col_mean, 0)
      }
    }
  }
  
  Z_samples_all <- array(NA, dim = c(n_chains, n_iter, N, K))
  W_samples_all <- array(NA, dim = c(n_chains, n_iter, P, K))
  mu_samples_all <- array(NA, dim = c(n_chains, n_iter, P))
  sigma2_samples_all <- matrix(NA, n_chains, n_iter)
  tau_g2_samples_all <- array(NA, dim = c(n_chains, n_iter, G))
  lambda2_samples_all <- matrix(NA, n_chains, n_iter)
  
  for (chain in 1:n_chains) {
    cat("Running chain", chain, "of", n_chains, "...\n")
    
    set.seed(123 + chain * 1000)
    
    mu <- colMeans(X_work, na.rm = TRUE)
    X_centered <- sweep(X_work, 2, mu)
    svd_init <- svd(X_centered, nu = min(N, K), nv = min(P, K))
    
    d_k <- svd_init$d[1:K]
    U_k <- svd_init$u[, 1:K, drop = FALSE]
    V_k <- svd_init$v[, 1:K, drop = FALSE]
    
    sqrt_d <- diag(sqrt(pmax(d_k / N, 1e-8)), nrow = K, ncol = K)
    Z <- U_k %*% sqrt_d * sqrt(N) + matrix(rnorm(N * K, 0, 0.1), N, K)
    W <- V_k %*% sqrt_d + matrix(rnorm(P * K, 0, 0.1), P, K)
    
    sigma2 <- max(var(as.vector(X_work - Z %*% t(W)), na.rm = TRUE), 0.01)
    tau_g2 <- rep(1, G)
    lambda2 <- 1
    
    if (!is.na(fix_sigma2)) {
      sigma2 <- fix_sigma2
    }
    if (isTRUE(no_shrinkage)) {
      tau_g2 <- rep(1e6, G)
      lambda2 <- 0
    }
    
    for (iter in 1:n_iter) {
      if (nrow(missing_idx) > 0) {
        X_pred <- Z %*% t(W) + matrix(rep(mu, each = N), nrow = N)
        X_work[missing_idx] <- rnorm(nrow(missing_idx), mean = X_pred[missing_idx], sd = sqrt(sigma2))
      }
      
      if (use_block_sampling && iter %% 3 == 1 && K > 1) {
        block_result <- update_W_Z_blocks(X_work, W, Z, mu, sigma2, group_id, G, K, no_shrinkage)
        W <- block_result$W
        Z <- block_result$Z
      } else {
        Z <- sample_z(W, X_work, mu, sigma2, N, K)
        X_residual_W <- X_work - matrix(rep(mu, each = N), nrow = N)
        W <- sample_group_lasso_loadings(X_residual_W, Z, sigma2, tau_g2, group_id, G)
      }
      
      for (k in 1:K) {
        w_sum <- sum(W[, k], na.rm = TRUE)
        if (is.finite(w_sum) && w_sum < 0) {
          W[, k] <- -W[, k]
          Z[, k] <- -Z[, k]
        }
      }
      
      prior_precision_mu <- 1 / mu_prior_var
      data_precision_mu <- N / sigma2
      post_var_mu <- 1 / (data_precision_mu + prior_precision_mu)
      for (j in 1:P) {
        residual_sum_mu <- sum(X_work[, j] - Z %*% W[j, ])
        post_mean_mu <- post_var_mu * (residual_sum_mu / sigma2)
        mu[j] <- rnorm(1, mean = post_mean_mu, sd = sqrt(post_var_mu))
      }
      
      if (is.na(fix_sigma2)) {
        X_pred_sigma <- Z %*% t(W) + matrix(rep(mu, each = N), nrow = N)
        sse <- sum((X_work - X_pred_sigma)^2)
        
        prior_penalty_sigma <- 0
        total_params_w <- 0
        for (g in 1:G) {
          vars_in_group <- which(group_id == g)
          if (length(vars_in_group) > 0) {
            W_g <- W[vars_in_group, , drop = FALSE]
            prior_penalty_sigma <- prior_penalty_sigma + sum(W_g^2) / tau_g2[g]
            total_params_w <- total_params_w + length(W_g)
          }
        }
        shape_post_sigma <- sigma_a + (N * P + total_params_w) / 2
        scale_post_sigma <- sigma_b + (sse + prior_penalty_sigma) / 2
        sigma2 <- 1 / rgamma(1, shape = shape_post_sigma, rate = scale_post_sigma)
        if (!is.finite(sigma2) || sigma2 < 1e-9) sigma2 <- 1.0
      } else {
        sigma2 <- fix_sigma2
      }
      
      if (!isTRUE(no_shrinkage)) {
        tau_g2 <- sample_group_scales_gig(W, group_id, sigma2, lambda2, G, K)
        lambda2 <- sample_global_shrinkage(tau_g2, group_sizes_vec, lambda_a, lambda_b, K)
      }
      
      if (iter > burn_in) {
        idx <- iter - burn_in
        Z_samples_all[chain, idx, , ] <- Z
        W_samples_all[chain, idx, , ] <- W
        mu_samples_all[chain, idx, ] <- mu
        sigma2_samples_all[chain, idx] <- sigma2
        tau_g2_samples_all[chain, idx, ] <- tau_g2
        lambda2_samples_all[chain, idx] <- lambda2
      }
    }
  }
  
  post_burn_iters <- n_iter - burn_in
  
  aligned_results <- align_factor_chains(
    W_samples_all[, 1:post_burn_iters, , ],
    Z_samples_all[, 1:post_burn_iters, , ]
  )
  W_samples_aligned <- aligned_results$W_aligned
  Z_samples_aligned <- aligned_results$Z_aligned
  
  Z_mean <- apply(Z_samples_aligned, c(3, 4), mean)
  W_mean <- apply(W_samples_aligned, c(3, 4), mean)
  mu_mean <- apply(mu_samples_all[, 1:post_burn_iters, ], 2, mean)
  sigma2_mean <- mean(sigma2_samples_all[, 1:post_burn_iters])
  
  convergence <- calculate_convergence_diagnostics(
    sigma2_samples_all[, 1:post_burn_iters, drop = FALSE], 
    tau_g2_samples_all[, 1:post_burn_iters, , drop = FALSE], 
    lambda2_samples_all[, 1:post_burn_iters, drop = FALSE],
    W_samples_aligned, Z_samples_aligned
  )
  
  return(list(
    Z_mean = Z_mean, W_mean = W_mean, mu_mean = mu_mean, sigma2_mean = sigma2_mean,
    convergence = convergence,
    W_samples_aligned = W_samples_aligned, Z_samples_aligned = Z_samples_aligned,
    mu_samples = mu_samples_all[, 1:post_burn_iters, , drop = FALSE],
    sigma2_samples = sigma2_samples_all[, 1:post_burn_iters, drop = FALSE],
    tau_g2_samples = tau_g2_samples_all[, 1:post_burn_iters, , drop = FALSE],
    lambda2_samples = lambda2_samples_all[, 1:post_burn_iters, drop = FALSE]
  ))
}

align_factor_chains <- function(W_chains, Z_chains, method = "procrustes") {
  # Ensure W_chains and Z_chains are 4-dimensional arrays
  if (length(dim(W_chains)) == 3) {
    # If K=1, add the K dimension back
    W_chains <- array(W_chains, dim = c(dim(W_chains), 1))
  }
  if (length(dim(Z_chains)) == 3) {
    # If K=1, add the K dimension back
    Z_chains <- array(Z_chains, dim = c(dim(Z_chains), 1))
  }
  
  n_chains <- dim(W_chains)[1]
  n_iter <- dim(W_chains)[2]
  K <- dim(W_chains)[4]
  
  if (n_chains < 2) {
    return(list(W_aligned = W_chains, Z_aligned = Z_chains))
  }
  
  W_aligned <- W_chains
  Z_aligned <- Z_chains
  
  W_means <- apply(W_chains, c(1, 3, 4), mean)
  
  for (chain in 2:n_chains) {
    W_ref <- W_means[1, , ]
    W_curr_mean <- W_means[chain, , ]
    
    if (K > 1) {
      cross_cov <- t(W_ref) %*% W_curr_mean
      if (rcond(cross_cov) > 1e-12 && !any(is.na(cross_cov))) {
        svd_result <- svd(cross_cov)
        R_opt <- svd_result$v %*% t(svd_result$u)
        
        if (det(R_opt) < 0) {
          svd_result$v[, K] <- -svd_result$v[, K]
          R_opt <- svd_result$v %*% t(svd_result$u)
        }
      } else {
        R_opt <- diag(K)
      }
    } else {
      sign_corr <- sign(sum(W_ref * W_curr_mean))
      if(is.na(sign_corr) || sign_corr == 0) sign_corr <- 1
      R_opt <- matrix(sign_corr)
    }
    
    for (iter in 1:n_iter) {
      W_aligned[chain, iter, , ] <- W_chains[chain, iter, , ] %*% R_opt
      Z_aligned[chain, iter, , ] <- Z_chains[chain, iter, , ] %*% R_opt
    }
  }
  
  return(list(W_aligned = W_aligned, Z_aligned = Z_aligned))
}

calculate_convergence_diagnostics <- function(
  sigma2_chains,
  tau_g_chains,
  lambda2_chains,
  W_chains,
  Z_chains,
  save_trace_plots = FALSE,
  trace_pdf_path = NULL,
  title_prefix = NULL
) {

  calculate_rhat <- function(chains_matrix) {
    n_chains <- nrow(chains_matrix)
    n_iter <- ncol(chains_matrix)
    if (n_iter < 20 || n_chains < 2) return(NA_real_)

    split_chains <- function(x) {
      half <- floor(ncol(x) / 2)
      if (half < 10) return(x)
      rbind(x[, 1:half, drop = FALSE], x[, (half + 1):ncol(x), drop = FALSE])
    }

    compute_rhat_core <- function(x) {
      m <- nrow(x); n <- ncol(x)
      chain_means <- rowMeans(x)
      chain_vars <- apply(x, 1, var)
      W <- mean(chain_vars)
      B <- n * var(chain_means)
      var_plus <- ((n - 1) / n) * W + (1 / n) * B
      if (!is.finite(W) || W <= 1e-12) return(Inf)
      sqrt(var_plus / W)
    }

    x <- split_chains(chains_matrix)
    rhat_basic <- compute_rhat_core(x)

    # Rank-normalized R-hat
    m <- nrow(x); n <- ncol(x)
    pooled <- as.vector(t(x))
    finite_mask <- is.finite(pooled)
    ranks <- rank(pooled[finite_mask], ties.method = "average")
    z <- rep(NA_real_, length(pooled))
    z[finite_mask] <- qnorm((ranks - 0.5) / sum(finite_mask))
    z_mat <- t(matrix(z, nrow = n, ncol = m))
    rhat_rank <- compute_rhat_core(z_mat)

    # Folded rank-normalized R-hat
    med <- median(pooled[finite_mask])
    folded <- abs(pooled - med)
    f_ranks <- rank(folded[finite_mask], ties.method = "average")
    fz <- rep(NA_real_, length(pooled))
    fz[finite_mask] <- qnorm((f_ranks - 0.5) / sum(finite_mask))
    fz_mat <- t(matrix(fz, nrow = n, ncol = m))
    rhat_folded <- compute_rhat_core(fz_mat)

    rhat <- max(rhat_basic, rhat_rank, rhat_folded, na.rm = TRUE)
    if (!is.finite(rhat)) rhat <- 2.0
    rhat
  }
  
  print_rhat <- function(name, value) {
    cat(sprintf("   %-30s R-hat: %.4f\n", name, value))
  }

  cat("=== CONVERGENCE ASSESSMENT ===\n\n")
  all_rhats <- c()
  
  cat("1. MODEL-WIDE PARAMETERS:\n")
  sigma2_rhat <- calculate_rhat(sigma2_chains)
  print_rhat("Noise Variance (sigma2)", sigma2_rhat)
  lambda2_rhat <- calculate_rhat(lambda2_chains)
  print_rhat("Global Shrinkage (lambda2)", lambda2_rhat)
  all_rhats <- c(sigma2_rhat, lambda2_rhat)
  
  cat("\n2. GROUP-SPECIFIC SCALES (tau_g):\n")
  num_tau_test <- min(5, dim(tau_g_chains)[3])
  tau_rhats <- sapply(1:num_tau_test, function(g) {
    rhat <- calculate_rhat(tau_g_chains[, , g])
    print_rhat(paste0("Group Scale (tau_", g, ")"), rhat)
    return(rhat)
  })
  all_rhats <- c(all_rhats, tau_rhats)
  
  cat("\n3. FACTOR LOADINGS (W):\n")
  P <- dim(W_chains)[3]
  K <- dim(W_chains)[4]
  w_indices <- list(c(1,1), c(ceiling(P/2), 1), c(P,1), c(ceiling(P/4), K), c(P, K))
  W_rhats <- sapply(w_indices, function(idx) {
    if (idx[1] <= P && idx[2] <= K) {
      rhat <- calculate_rhat(W_chains[, , idx[1], idx[2]])
      print_rhat(paste0("Loading W[", idx[1], ",", idx[2], "]"), rhat)
      return(rhat)
    }
    return(NA)
  })
  all_rhats <- c(all_rhats, W_rhats)
  
  cat("\n4. LATENT FACTORS (Z):\n")
  N <- dim(Z_chains)[3]
  z_indices <- list(c(1,1), c(ceiling(N/2), 1), c(N,1), c(ceiling(N/4), K), c(N, K))
  Z_rhats <- sapply(z_indices, function(idx) {
    if (idx[1] <= N && idx[2] <= K) {
      rhat <- calculate_rhat(Z_chains[, , idx[1], idx[2]])
      print_rhat(paste0("Factor Z[", idx[1], ",", idx[2], "]"), rhat)
      return(rhat)
    }
    return(NA)
  })
  all_rhats <- c(all_rhats, Z_rhats)
  
  all_rhats <- all_rhats[!is.na(all_rhats)]
  max_rhat <- max(all_rhats, na.rm = TRUE)
  prop_converged_excellent <- mean(all_rhats < 1.1, na.rm = TRUE)
  
  cat("\n=== CONVERGENCE SUMMARY ===\n")
  cat(sprintf("   Total parameters tested: %d\n", length(all_rhats)))
  cat(sprintf("   Max R-hat: %.4f\n", max_rhat))
  cat(sprintf("   Convergence (< 1.1): %.1f%%\n", 100 * prop_converged_excellent))
  
  if (max_rhat < 1.1) {
    cat("\nRECOMMENDATION: EXCELLENT. All tested parameters have converged.\n")
  } else if (max_rhat < 1.2) {
    cat("\nRECOMMENDATION: GOOD. Most parameters have converged. Consider more iterations.\n")
  } else {
    cat("\nRECOMMENDATION: POOR. Convergence issues detected. Increase iterations or check model.\n")
  }
  
  if (save_trace_plots && !is.null(trace_pdf_path)) {
    to_mcmc_list <- function(chains_matrix) {
      coda::mcmc.list(lapply(seq_len(nrow(chains_matrix)), function(ch) coda::mcmc(chains_matrix[ch, ])))
    }

    dir.create(dirname(trace_pdf_path), recursive = TRUE, showWarnings = FALSE)
    # Taller page for better readability; 3:4 aspect
    pdf(trace_pdf_path, width = 12, height = 16)
    on.exit(try({ dev.off() }, silent = TRUE), add = TRUE)

    plots <- list()
    plots[[length(plots) + 1]] <- list(obj = to_mcmc_list(sigma2_chains), title = "sigma2")
    plots[[length(plots) + 1]] <- list(obj = to_mcmc_list(lambda2_chains), title = "lambda2")

    num_tau_test <- min(5, dim(tau_g_chains)[3])
    for (g in seq_len(num_tau_test)) {
      plots[[length(plots) + 1]] <- list(
        obj = to_mcmc_list(tau_g_chains[, , g]),
        title = paste0("tau_g[", g, "]")
      )
    }

    P <- dim(W_chains)[3]
    K <- dim(W_chains)[4]
    w_indices <- list(c(1, 1), c(ceiling(P / 2), 1), c(P, 1), c(ceiling(P / 4), K), c(P, K))
    for (idx in w_indices) {
      if (idx[1] <= P && idx[2] <= K) {
        plots[[length(plots) + 1]] <- list(
          obj = to_mcmc_list(W_chains[, , idx[1], idx[2]]),
          title = paste0("W[", idx[1], ",", idx[2], "]")
        )
      }
    }

    N <- dim(Z_chains)[3]
    z_indices <- list(c(1, 1), c(ceiling(N / 2), 1), c(N, 1), c(ceiling(N / 4), K), c(N, K))
    for (idx in z_indices) {
      if (idx[1] <= N && idx[2] <= K) {
        plots[[length(plots) + 1]] <- list(
          obj = to_mcmc_list(Z_chains[, , idx[1], idx[2]]),
        title = paste0("Z[", idx[1], ",", idx[2], "]")
        )
      }
    }

    # Use a fixed 3x2 grid per page; additional plots overflow to new pages automatically
    n_cols <- 2
    n_rows <- 3
    op <- par(mfrow = c(n_rows, n_cols), mar = c(3.5, 3.5, 2.2, 1), cex.lab = 1.1, cex.axis = 0.95)
    on.exit(try({ par(op) }, silent = TRUE), add = TRUE)

    for (pl in plots) {
      coda::traceplot(pl$obj, main = pl$title)
    }
  }

  return(invisible(list(
    max_rhat = max_rhat,
    all_rhats = all_rhats,
    prop_converged = prop_converged_excellent,
    trace_pdf = if (save_trace_plots) trace_pdf_path else NULL
  )))
}

evaluate_factor_recovery <- function(W_true, W_est) {
  pro <- vegan::procrustes(X = W_true, Y = W_est, symmetric = TRUE)
  W_est_aligned <- W_est %*% pro$rotation
  
  cor_loadings <- cor(as.vector(W_true), as.vector(W_est_aligned))
  rmse_loadings <- sqrt(mean((W_true - W_est_aligned)^2))
  
  return(list(
    W_aligned = W_est_aligned,
    correlation = cor_loadings,
    rmse = rmse_loadings,
    procrustes_ss = pro$ss
  ))
}

plot_W_comparison <- function(W_true, W_est_aligned, title_suffix = "", save_plots = TRUE, output_dir = "bpca_missing_lasso_plots") {
  cat("=== CREATING W MATRIX HEATMAP COMPARISON ===\n\n")
  
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  cat("1. Normalization Tests:\n")
  
  norm_true <- norm(W_true, type = "F")
  norm_est <- norm(W_est_aligned, type = "F")
  norm_ratio <- norm_est / norm_true
  
  cat("   Frobenius norm (true):", format(round(norm_true, 4), nsmall = 4), "\n")
  cat("   Frobenius norm (estimated):", format(round(norm_est, 4), nsmall = 4), "\n")
  cat("   Norm ratio (est/true):", format(round(norm_ratio, 4), nsmall = 4), "\n")
  
  col_norms_true <- apply(W_true, 2, function(x) norm(x, type = "2"))
  col_norms_est <- apply(W_est_aligned, 2, function(x) norm(x, type = "2"))
  
  cat("   Column norms (true):", paste(format(round(col_norms_true, 4), nsmall = 4), collapse = ", "), "\n")
  cat("   Column norms (estimated):", paste(format(round(col_norms_est, 4), nsmall = 4), collapse = ", "), "\n")

  assess_norm <- function(r) {
    if (!is.finite(r)) return("undetermined")
    if (r > 0.9 && r < 1.1) return("excellent recovery (scale matched)")
    if (r > 0.75 && r < 1.25) return("good recovery (minor scale drift)")
    return("poor recovery (scale mismatch)")
  }
  norm_assessment <- assess_norm(norm_ratio)
  cat("   Interpretation (norm ratio):", norm_assessment, "\n\n")
  
  cat("2. Creating Heatmaps...\n")
  
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
  
  W_diff <- W_est_aligned - W_true
  rownames(W_diff) <- matrix_rownames
  colnames(W_diff) <- matrix_colnames
  
  all_values <- c(as.vector(W_true_labeled), as.vector(W_est_labeled))
  min_val <- min(all_values)
  max_val <- max(all_values)
  color_breaks <- seq(min_val, max_val, length.out = 100)
  my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
  
  diff_max <- max(abs(W_diff))
  diff_breaks <- seq(-diff_max, diff_max, length.out = 100)
  diff_colors <- colorRampPalette(c("darkred", "white", "darkblue"))(99)
  
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
  
  if (save_plots) {
    method_name <- gsub("[^A-Za-z0-9]", "_", title_suffix)
    if (method_name == "") method_name <- "standard"
    
    pdf(file.path(output_dir, paste0("W_comparison_combined", method_name, ".pdf")), width = 24, height = 10)
    grid.arrange(p_true$gtable, p_est$gtable, p_diff$gtable, ncol = 3)
    dev.off()
    
    cat("  Saved W comparison plots to:", output_dir, "\n")
  } else {
    grid.arrange(p_true$gtable, p_est$gtable, p_diff$gtable, ncol = 3)
  }
  
  return(list(
    norm_ratio = norm_ratio,
    col_norms_ratio = col_norms_est / col_norms_true,
    max_abs_diff = max(abs(W_diff)),
    norm_assessment = norm_assessment,
    plots = list(true = p_true, estimated = p_est, difference = p_diff)
  ))
}

debug_scale_correction <- function(W_chains, Z_chains, W_true = NULL) {
  cat("=== DEBUGGING SCALE CORRECTION ===\n\n")
  W_est_original <- apply(W_chains, c(3, 4), mean)
  Z_est <- apply(Z_chains, c(3, 4), mean)

  cat("1. POSTERIOR MEANS DIMENSIONS:\n")
  cat(sprintf("   W_est: %d x %d\n", nrow(W_est_original), ncol(W_est_original)))
  cat(sprintf("   Z_est: %d x %d\n", nrow(Z_est), ncol(Z_est)))

  cat("\n2. Z COLUMNS STATISTICS:\n")
  K <- ncol(Z_est)
  z_stats <- data.frame(
    Factor = 1:K,
    Mean = numeric(K),
    SD = numeric(K),
    Min = numeric(K),
    Max = numeric(K),
    N_finite = numeric(K),
    Will_scale = logical(K)
  )
  for (k in seq_len(K)) {
    z_col <- Z_est[, k]
    z_stats$Mean[k] <- mean(z_col, na.rm = TRUE)
    z_stats$SD[k] <- sd(z_col, na.rm = TRUE)
    z_stats$Min[k] <- min(z_col, na.rm = TRUE)
    z_stats$Max[k] <- max(z_col, na.rm = TRUE)
    z_stats$N_finite[k] <- sum(is.finite(z_col))
    z_stats$Will_scale[k] <- is.finite(z_stats$SD[k]) && z_stats$SD[k] > 1e-6
    cat(sprintf(
      "   Factor %d: mean=%.4f, sd=%.4f, range=[%.4f, %.4f], finite=%d/%d, will_scale=%s\n",
      k, z_stats$Mean[k], z_stats$SD[k], z_stats$Min[k], z_stats$Max[k],
      z_stats$N_finite[k], nrow(Z_est), z_stats$Will_scale[k]
    ))
  }

  cat("\n3. W COLUMNS INITIAL STATE:\n")
  w_norms_before <- apply(W_est_original, 2, function(x) norm(x, "2"))
  cat("   W column norms (before):", paste(round(w_norms_before, 4), collapse = ", "), "\n")

  cat("\n4. APPLYING SCALE CORRECTION:\n")
  W_est_corrected <- W_est_original
  scale_factors_applied <- numeric(K)
  for (k in seq_len(K)) {
    s <- sd(Z_est[, k], na.rm = TRUE)
    cat(sprintf("   Factor %d: Z sd = %.6f", k, s))
    if (is.finite(s) && s > 1e-6) {
      Z_est[, k] <- Z_est[, k] / s
      W_est_corrected[, k] <- W_est_corrected[, k] * s
      scale_factors_applied[k] <- s
      cat(sprintf(" -> APPLIED (scale factor = %.4f)\n", s))
    } else {
      scale_factors_applied[k] <- 1.0
      cat(" -> SKIPPED (invalid sd)\n")
    }
  }

  cat("\n5. BEFORE/AFTER COMPARISON:\n")
  w_norms_after <- apply(W_est_corrected, 2, function(x) norm(x, "2"))
  z_sds_after <- apply(Z_est, 2, sd, na.rm = TRUE)
  cat("   W column norms before:", paste(round(w_norms_before, 4), collapse = ", "), "\n")
  cat("   W column norms after: ", paste(round(w_norms_after, 4), collapse = ", "), "\n")
  cat("   Z column sds after:   ", paste(round(z_sds_after, 4), collapse = ", "), "\n")
  cat("   Scale factors applied: ", paste(round(scale_factors_applied, 4), collapse = ", "), "\n")

  if (!is.null(W_true)) {
    norm_before <- norm(W_est_original, "F") / norm(W_true, "F")
    norm_after <- norm(W_est_corrected, "F") / norm(W_true, "F")
    cat("\n6. NORM RATIO ANALYSIS:\n")
    cat(sprintf("   True W Frobenius norm: %.4f\n", norm(W_true, "F")))
    cat(sprintf("   Norm ratio before correction: %.4f\n", norm_before))
    cat(sprintf("   Norm ratio after correction:  %.4f\n", norm_after))
    cat(sprintf("   Change: %.4f -> %.4f (%+.4f)\n", norm_before, norm_after, norm_after - norm_before))
    if (norm_after < norm_before) {
      cat("   WARNING: Norm ratio got WORSE after correction!\n")
    }
  }

  cat("\n7. NUMERICAL HEALTH CHECK:\n")
  w_has_inf <- any(!is.finite(W_est_corrected))
  z_has_inf <- any(!is.finite(Z_est))
  w_max_abs <- max(abs(W_est_corrected), na.rm = TRUE)
  z_max_abs <- max(abs(Z_est), na.rm = TRUE)
  cat(sprintf("   W has non-finite values: %s\n", w_has_inf))
  cat(sprintf("   Z has non-finite values: %s\n", z_has_inf))
  cat(sprintf("   W max absolute value: %.4f\n", w_max_abs))
  cat(sprintf("   Z max absolute value: %.4f\n", z_max_abs))
  if (w_max_abs > 1e3 || z_max_abs > 1e3) {
    cat("   WARNING: Very large values detected!\n")
  }

  return(invisible(list(
    W_original = W_est_original,
    W_corrected = W_est_corrected,
    Z_corrected = Z_est,
    scale_factors = scale_factors_applied,
    z_statistics = z_stats
  )))
}

debug_and_fix_scale_correction <- function(W_chains, Z_chains, W_true = NULL) {
  dbg <- debug_scale_correction(W_chains, Z_chains, W_true)
  if (any(dbg$scale_factors < 1e-6) || any(!is.finite(dbg$scale_factors))) {
    cat("\n=== DETECTED SCALE CORRECTION ISSUES, TRYING ALTERNATIVE ===\n")
    W_est <- apply(W_chains, c(3, 4), mean)
    for (k in seq_len(ncol(W_est))) {
      w_norm <- norm(W_est[, k], "2")
      if (is.finite(w_norm) && w_norm > 1e-6) {
        W_est[, k] <- W_est[, k] / w_norm
      }
    }
    cat("Applied alternative normalization based on W column norms\n")
    return(W_est)
  } else {
    cat("\nScale correction appears successful\n")
    return(dbg$W_corrected)
  }
}

diagnosis_report <- function(W_chains, Z_chains, sigma2_chains, tau_g_chains, lambda2_chains,
                             W_true = NULL, title_suffix = "", save_plots = TRUE,
                             output_dir = "bpca_missing_lasso_plots",
                             save_trace_plots = FALSE,
                             trace_pdf_path = NULL) {
  # 1) Convergence diagnostics (uses a subset across parameters for speed/readability)
  convergence <- calculate_convergence_diagnostics(
    sigma2_chains = sigma2_chains,
    tau_g_chains = tau_g_chains,
    lambda2_chains = lambda2_chains,
    W_chains = W_chains,
    Z_chains = Z_chains,
    save_trace_plots = save_trace_plots,
    trace_pdf_path = if (!is.null(trace_pdf_path)) trace_pdf_path else file.path(output_dir, "comprehensive_mcmc_traces.pdf"),
    title_prefix = NULL
  )

  # 2) Factor recovery
  recovery <- NULL
  plot_info <- NULL
  if (!is.null(W_true)) {
    # ======== DISABLE SCALE CORRECTION FOR DEBUGGING (KEEP CALL COMMENTED) ========
    # W_est <- debug_and_fix_scale_correction(W_chains, Z_chains, W_true)
    # Using raw posterior mean of W for evaluation/plots to avoid masking scale issues
    W_est <- apply(W_chains, c(3, 4), mean)
    # ======== END DISABLE BLOCK ========
    # Posterior mean of W across chains and iterations
    # W_est <- apply(W_chains, c(3, 4), mean)
    recovery <- evaluate_factor_recovery(W_true = W_true, W_est = W_est)

    # 3) Plot comparisons
    plot_info <- plot_W_comparison(
      W_true = W_true,
      W_est_aligned = recovery$W_aligned,
      title_suffix = title_suffix,
      save_plots = save_plots,
      output_dir = output_dir
    )
  }

  return(invisible(list(
    convergence = convergence,
    recovery = recovery,
    plots = plot_info
  )))
}