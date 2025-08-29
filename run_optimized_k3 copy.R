source("bpca_missing_lasso_concise copy.R")
source("bayesian_group_lasso_true copy.R")

set.seed(123)

# Output directory for this low-dimensional setup
output_dir <- "p30k3_fix_w11"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Always log this run
log_path <- file.path(output_dir, "run.log")
cat(sprintf("\n===== Run started: %s =====\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    file = log_path, append = TRUE)
con_out <- file(log_path, open = "a")
sink(con_out, split = TRUE)
sink(con_out, type = "message")
on.exit({
  try(sink(NULL), silent = TRUE)
  try(sink(NULL, type = "message"), silent = TRUE)
  try(close(con_out), silent = TRUE)
  cat(sprintf("===== Run ended: %s =====\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file = log_path, append = TRUE)
}, add = TRUE)

sim_data <- simulate_data_with_groups(
  n = 200,
  p = 30,
  G = 30,
  K = 3,
  # One variable per group (G = p). Activate only a subset to leave blank groups (sparser).
  active_map = list(
    "1" = 1:5,
    "2" = 6:10,
    "3" = 11:15
  ),
  sigma = 0.2
)

X_data_complete <- sim_data$X_data
W_true <- sim_data$W_true

set.seed(123)
missing_rate <- 0.10
n_missing <- round(nrow(X_data_complete) * ncol(X_data_complete) * missing_rate)
missing_indices <- sample(seq_along(X_data_complete), n_missing)
X_data_missing <- X_data_complete
X_data_missing[missing_indices] <- NA

result_lasso <- mcmc_bpca_group_lasso(
  X = X_data_missing,
  K = 3,
  G = 30,
  n_iter = 30000,
  burn_in = 15000,
  n_chains = 3,
  lambda_a = 0.1,
  lambda_b = 0.1,
  # ======== ISOLATION TOGGLES (ADDED) ========
  fix_sigma2 = 0.04,
  no_shrinkage = TRUE
)

# Convergence diagnostics and save
diag <- diagnosis_report(
  W_chains = result_lasso$W_samples_aligned,
  Z_chains = result_lasso$Z_samples_aligned,
  sigma2_chains = result_lasso$sigma2_samples,
  tau_g_chains = result_lasso$tau_g2_samples,
  lambda2_chains = result_lasso$lambda2_samples,
  W_true = W_true,
  save_plots = TRUE,
  output_dir = output_dir,
  save_trace_plots = TRUE,
  trace_pdf_path = file.path(output_dir, "comprehensive_mcmc_traces.pdf")
)

conv <- diag$convergence

saveRDS(result_lasso, file.path(output_dir, "mcmc_result.rds"))
saveRDS(conv, file.path(output_dir, "convergence_summary.rds"))

if (!is.null(conv$max_rhat)) {
  recovery <- diag$recovery
  plot_info <- diag$plots
  conv_df <- data.frame(
    max_rhat = round(conv$max_rhat, 4),
    prop_converged_lt_1_1 = round(conv$prop_converged, 4),
    w_recovery_correlation = if (!is.null(recovery)) round(recovery$correlation, 4) else NA_real_,
    w_recovery_rmse = if (!is.null(recovery)) round(recovery$rmse, 4) else NA_real_,
    w_recovery_procrustes_ss = if (!is.null(recovery)) round(recovery$procrustes_ss, 4) else NA_real_,
    w_norm_ratio = if (!is.null(plot_info)) round(plot_info$norm_ratio, 4) else NA_real_,
    w_norm_ratio_interpretation = if (!is.null(plot_info)) plot_info$norm_assessment else NA_character_,
    w_max_abs_diff = if (!is.null(plot_info)) round(plot_info$max_abs_diff, 4) else NA_real_,
    w_col_norms_ratio_mean = if (!is.null(plot_info)) round(mean(plot_info$col_norms_ratio), 4) else NA_real_
  )
  write.csv(conv_df, file.path(output_dir, "convergence_summary.csv"), row.names = FALSE)

  log_path <- file.path(output_dir, "run.log")
  log_line <- sprintf(
    "%s | max_rhat=%.3f, prop_converged=%.3f, corr=%.3f, rmse=%.4f, norm_ratio=%.3f, max_abs_diff=%.3f\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    conv$max_rhat,
    conv$prop_converged,
    if (!is.null(recovery)) recovery$correlation else NA_real_,
    if (!is.null(recovery)) recovery$rmse else NA_real_,
    if (!is.null(plot_info)) plot_info$norm_ratio else NA_real_,
    if (!is.null(plot_info)) plot_info$max_abs_diff else NA_real_
  )
  cat(log_line, file = log_path, append = TRUE)
}