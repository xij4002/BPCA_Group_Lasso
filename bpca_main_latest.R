library(mvtnorm)
library(GIGrvg)
library(MCMCpack)

sample_gig <- function(p, a, b) {
  if (!is.finite(p) || !is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    return(1.0) 
  }
  
  tryCatch({
    result <- GIGrvg::rgig(1, lambda = p, chi = b, psi = a)
    if (!is.finite(result) || result <= 1e-9 || result > 1e6) {
      return(1.0)
    }
    return(result)
  }, error = function(e) {
    return(1.0)
  })
}

sample_group_lasso_loadings <- function(X_residual, Z, sigma2, tau_g_squared, group_id, G) {
  P <- ncol(X_residual)
  K <- ncol(Z)
  W <- matrix(0, nrow = P, ncol = K)
  Z_T_Z <- crossprod(Z)
  
  for (g in 1:G) {
    vars_in_group <- which(group_id == g)
    d_g <- length(vars_in_group)
    if (d_g == 0) next
    X_residual_g <- X_residual[, vars_in_group, drop = FALSE]
    
    Omega <- Z_T_Z + (1 / tau_g_squared[g]) * diag(K)
    min_eig <- min(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values)
    if (min_eig < 1e-10) {
      Omega <- Omega + (1e-8 - min_eig) * diag(K)
    }
    B_g <- crossprod(X_residual_g, Z)
  
    tryCatch({
      L_omega <- chol(Omega)  # Omega = L_omega^T L_omega
      # Omega * mu^T = B_g^T
      mu_W_g <- t(backsolve(L_omega, backsolve(L_omega, t(B_g), transpose = TRUE)))
      
      # W_g ~ MN(mu_W_g, sigma2 * I, Omega^{-1})
      Z_std <- matrix(rnorm(d_g * K), d_g, K)
      
      # W_g = mu_W_g + sqrt(sigma2) * Z_std * L_omega^{-1}
      L_omega_inv <- backsolve(L_omega, diag(K))
      W_g_new <- mu_W_g + sqrt(sigma2) * Z_std %*% t(L_omega_inv)
      
      W[vars_in_group, ] <- W_g_new
      
    }, error = function(e) {
      warning(paste("Group", g, ": Using regularized fallback"))
      A_g <- Omega + 1e-6 * diag(K)
      A_g_inv <- solve(A_g)
      mu_W_g <- B_g %*% A_g_inv
      Sigma_W_g <- sigma2 * kronecker(A_g_inv, diag(d_g))
      vec_mu_W_g <- as.vector(mu_W_g)
      vec_W_g_new <- mvtnorm::rmvnorm(1, mean = vec_mu_W_g, sigma = Sigma_W_g)
      W_g_new <- matrix(vec_W_g_new, nrow = d_g, ncol = K, byrow = FALSE)
      W[vars_in_group, ] <- W_g_new
    })
  }
  return(W)
}

sample_group_scales_gig <- function(W, group_id, sigma2, lambda2, G, K) {
  tau_g_squared <- numeric(G)
  
  for (g in 1:G) {
    vars_in_group <- which(group_id == g)
    if (length(vars_in_group) == 0) {
      tau_g_squared[g] <- 1.0
      next
    }
    
    W_g <- W[vars_in_group, , drop = FALSE]
  
    p_gig <- 1/2
    a_gig <- lambda2
    b_gig <- sum(W_g^2) / sigma2
    
    tau_g_squared[g] <- sample_gig(p = p_gig, a = a_gig, b = b_gig)
  }
  
  return(tau_g_squared)
}

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
