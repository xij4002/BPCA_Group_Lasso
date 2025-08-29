# Load function libraries
source("bpca_missing_lasso_concise_fixw11_robust copy.R")
source("bayesian_group_lasso_true_fixw11_robust.R")
library(pheatmap)

# Generates a consistent, un-clustered heatmap for diagnostics
generate_diagnostic_heatmap <- function(W_samples, title, output_path) {
  n_chains <- dim(W_samples)[1]
  n_iter <- dim(W_samples)[2]
  P <- dim(W_samples)[3]
  K <- dim(W_samples)[4]

  W_flat <- array(W_samples, dim = c(n_chains * n_iter, P * K))
  
  cat(sprintf("Calculating correlation matrix for '%s'...\n", title))
  cor_matrix <- cor(W_flat)
  
  pdf(output_path, width = 12, height = 12)
  pheatmap::pheatmap(
    cor_matrix,
    main = title,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    gaps_row = c(P, P*2),
    gaps_col = c(P, P*2),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-1, 1, length.out = 101)
  )
  dev.off()
  cat(sprintf("'%s' heatmap generated.\n", title))
}


# Finds the most influential loading for each factor to use as an anchor
find_anchor_loadings <- function(W_samples) {
  if (length(dim(W_samples)) == 4 && dim(W_samples)[1] > 1) {
    W_mean <- apply(W_samples, c(3, 4), mean)
  } else {
    W_mean <- apply(W_samples[1, , , , drop = FALSE], c(3, 4), mean)
  }
  
  K <- ncol(W_mean)
  anchor_indices <- vector("list", K)
  names(anchor_indices) <- 1:K
  
  cat("--- Identifying Anchor Loadings ---\n")
  for (k in 1:K) {
    anchor_row <- which.max(abs(W_mean[, k]))
    anchor_indices[[k]] <- anchor_row
    cat(sprintf("Anchor for Factor %d: W[%d, %d] (Value: %.4f)\n", k, anchor_row, k, W_mean[anchor_row, k]))
  }
  cat("------------------------------------\n")
  return(anchor_indices)
}


# --- Main Workflow ---

set.seed(123)
output_dir <- "p30k3_final_auto_constraints"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Logging Setup ---
log_path <- file.path(output_dir, "run.log")
cat(sprintf("\n===== Auto-Constraint Run started: %s =====\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = log_path, append = TRUE)
con_out <- file(log_path, open = "a")
sink(con_out, split = TRUE)
sink(con_out, type = "message")
on.exit({
  try(sink(NULL), silent = TRUE)
  try(sink(NULL, type = "message"), silent = TRUE)
  try(close(con_out), silent = TRUE)
  cat(sprintf("===== Run ended: %s =====\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = log_path, append = TRUE)
}, add = TRUE)

# --- Data Simulation ---
sim_data <- simulate_data_with_groups(
  n = 200, p = 30, G = 30, K = 3,
  active_map = list("1" = 1:5, "2" = 6:10, "3" = 11:15),
  sigma = 0.2
)
X_data_complete <- sim_data$X_data
W_true <- sim_data$W_true
missing_rate <- 0.10
n_missing <- round(prod(dim(X_data_complete)) * missing_rate)
missing_indices <- sample(seq_along(X_data_complete), n_missing)
X_data_missing <- X_data_complete
X_data_missing[missing_indices] <- NA

# --- STEP 1: Exploratory MCMC Run ---
cat("\n--- STEP 1: Running Exploratory MCMC (1 chain, no constraints) ---\n")
exploratory_run <- mcmc_bpca_group_lasso(
  X = X_data_missing, K = 3, G = 30,
  n_iter = 60000, burn_in = 30000, n_chains = 1,
  fixed_loadings = NULL, no_shrinkage = TRUE, fix_sigma2 = 0.04, thinning = 10
)
cat("Exploratory run completed.\n")

generate_diagnostic_heatmap(
  W_samples = exploratory_run$W_samples_aligned,
  title = "Posterior Correlation (BEFORE Constraints)",
  output_path = file.path(output_dir, "posterior_correlation_BEFORE.pdf")
)

# --- STEP 2: Identify Anchors and Define Constraints ---
cat("\n--- STEP 2: Identifying Anchors and Constraint Values ---\n")
anchor_indices <- find_anchor_loadings(exploratory_run$W_samples_aligned)
W_exploratory_mean <- apply(exploratory_run$W_samples_aligned[1, , , ], c(2, 3), mean)
constraints_to_apply <- list()
for (k in 1:2) {  # Only fix factors 1 and 2
  anchor_row <- anchor_indices[[as.character(k)]]
  constraint_value <- W_exploratory_mean[anchor_row, k]
  constraints_to_apply[[k]] <- list(row = anchor_row, col = k, val = constraint_value)
  cat(sprintf("Constraint for final run: W[%d, %d] will be fixed to %.4f\n", anchor_row, k, constraint_value))
}
cat("Note: Factor 3 will be left unconstrained.\n")

# --- STEP 3: Final MCMC Run with Constraints ---
cat(sprintf("\n--- STEP 3: Running Final MCMC with %d constraints ---\n", length(constraints_to_apply)))
result_final <- mcmc_bpca_group_lasso(
  X = X_data_missing, K = 3, G = 30,
  n_iter = 100000, burn_in = 50000, n_chains = 3,
  fixed_loadings = constraints_to_apply,
  no_shrinkage = TRUE, fix_sigma2 = 0.04, thinning = 10
)
cat("Final constrained model MCMC run completed.\n")

generate_diagnostic_heatmap(
  W_samples = result_final$W_samples_aligned,
  title = "Posterior Correlation (AFTER Constraints)",
  output_path = file.path(output_dir, "posterior_correlation_AFTER.pdf")
)

# --- STEP 4: Final Diagnostics and Reporting ---
cat("\n--- STEP 4: Generating Final Diagnostic Report ---\n")
diag_report <- diagnosis_report(
  W_chains = result_final$W_samples_aligned,
  Z_chains = result_final$Z_samples_aligned,
  sigma2_chains = result_final$sigma2_samples,
  tau_g_chains = result_final$tau_g2_samples,
  lambda2_chains = result_final$lambda2_samples,
  W_true = W_true,
  title_suffix = " (Auto-Constraints)",
  save_plots = TRUE, output_dir = output_dir, save_trace_plots = TRUE,
  trace_pdf_path = file.path(output_dir, "mcmc_trace_plots_final.pdf")
)
cat("Full diagnostic report generated.\n")

saveRDS(result_final, file.path(output_dir, "mcmc_result_final.rds"))
saveRDS(diag_report, file.path(output_dir, "diagnostic_report_final.rds"))

if (!is.null(diag_report$convergence)) {
  conv <- diag_report$convergence
  recovery <- diag_report$recovery
  plot_info <- diag_report$plots
  
  summary_df <- data.frame(
    max_rhat = round(conv$max_rhat, 4),
    prop_converged_lt_1_1 = round(conv$prop_converged, 4),
    w_recovery_correlation = ifelse(!is.null(recovery), round(recovery$correlation, 4), NA),
    w_recovery_rmse = ifelse(!is.null(recovery), round(recovery$rmse, 4), NA),
    w_norm_ratio = ifelse(!is.null(plot_info), round(plot_info$norm_ratio, 4), NA)
  )
  write.csv(summary_df, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
}

cat("All processes completed successfully.\n")