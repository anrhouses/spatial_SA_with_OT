###############################################################################
# Sensitivity analysis with optimal transport

#1. Test several binning options accounting for the floor of 100 scenarios per bin
#2. SA-OT with the final bin and bootstrapping for CI
#3. Accompanying plots

#Author: Corinna Perchtold
###############################################################################
library(dplyr)
library(gsaot)
library(ggplot2)
library(patchwork)
library(future.apply) 

plan(multisession)          # use available cores; set workers = N to cap
SEED <- 20260624

load("01_Uncertainty_18k_preprocessed.R")

load("01_PCA_RMSE_log10_Cov95.R")

labelled <- labelled %>%
  filter(retain_metric %in% c("RMSE_log10","Coverage95"))


print(unique(labelled$retain_metric))
#shared configurations
THRESH_QUANTILE <- 0.95     # upper edge of the irrelevance zone
SOLVER          <- "sinkhorn"   
SINK_ITERS      <- 5e5      
SINK_MAXERR     <- 1e-4
BOOT_R          <- 200     # bootstrap replicates for Stage 2 CIs

# Stage 1: M-selection (fast, biased, not reported)
EPS_SELECT      <- 0.01
N_REPS_SELECT   <- 50
H_GRID          <- c(3, 5, 7, 9, 11, 13, 15,17,19,21,23,25)



sinkhorn_select <- list(epsilon = EPS_SELECT, numIterations = SINK_ITERS,
                        maxErr = SINK_MAXERR)


#check here uncertainty factors
uncertainty_factors <- labelled[, 2:6]

for (j in 1:3) uncertainty_factors[, j] <- as.factor(uncertainty_factors[, j])

BETA          <- as.matrix(pca_results_mean_rmse_log10_cov95$BETA)
n_pcs_to_test <- pca_results_mean_rmse_log10_cov95$n_pcs
N             <- nrow(uncertainty_factors)
n_fac         <- ncol(uncertainty_factors)

# Points-per-bin guard: gsaot recommends >= 100 points per partition set.
stopifnot(N / max(H_GRID) >= 100)
message(sprintf("Points per bin at largest H (%d): %.0f",
                max(H_GRID), N / max(H_GRID)))

threshold_quantile <- function(draws)
  as.numeric(quantile(draws, probs = THRESH_QUANTILE, names = FALSE))

# Parallel, reproducible dummy-draw threshold for one (y, M, solver) setting.
# Uniform dummy to match the uniform input draws.
dummy_threshold <- function(y, m_val, solver, solver_optns, n_reps, seed) {
  draws <- future_replicate(n_reps, {
    irrelevance_threshold(y = y, M = m_val, solver = solver,
                          solver_optns = solver_optns,
                          dummy_optns = list(distr = "runif"))$indices
  }, future.seed = seed)
  threshold_quantile(draws)
}


message("\n=== STAGE 1: H-selection sweep (", SOLVER, ") ===")

stage1 <- lapply(H_GRID, function(m_val) {
  message("  H = ", m_val)
  
  glob_idx <- ot_indices(x = uncertainty_factors, y = BETA, M = m_val,
                         solver = SOLVER, boot = FALSE,
                         solver_optns = sinkhorn_select)$indices
  
  set.seed(SEED)
  glob_thr <- dummy_threshold(BETA, m_val, SOLVER,
                              sinkhorn_select, N_REPS_SELECT, SEED)
  
  # Per-PC uses the exact 1d solver (no epsilon, always stable).
  pc <- lapply(1:n_pcs_to_test, function(i) {
    idx <- ot_indices_1d(x = uncertainty_factors, y = BETA[, i],
                         M = m_val, boot = FALSE)$indices
    thr <- dummy_threshold(BETA[, i], m_val, "1d", NULL,
                           N_REPS_SELECT, SEED + i)
    list(index = idx, threshold = thr)
  })
  
  data.frame(
    M      = m_val,
    Factor = rep(colnames(uncertainty_factors), n_pcs_to_test + 1),
    Index  = c(as.numeric(glob_idx),
               unlist(lapply(pc, `[[`, "index"))),
    Threshold = rep(c(glob_thr, unlist(lapply(pc, `[[`, "threshold"))),
                    each = n_fac),
    Type = rep(c("Global", paste0("PC", 1:n_pcs_to_test)), each = n_fac)
  )
}) %>% bind_rows()

band1 <- stage1 %>% group_by(M, Type) %>%
  summarise(Threshold = mean(Threshold), .groups = "drop")

plot<-ggplot() +
  geom_ribbon(data = band1, aes(M, ymin = 0, ymax = Threshold),
              fill = "firebrick", alpha = 0.15) +
  geom_line(data = band1, aes(M, Threshold),
            color = "firebrick4", linetype = "dotted", linewidth = 0.8) +
  geom_line(data = stage1, aes(M, Index, color = Factor, group = Factor),
            linewidth = 1) +
  geom_point(data = stage1, aes(M, Index, color = Factor), size = 2.5) +
  facet_wrap(~Type) +coord_cartesian(ylim = c(0, 1))+
  labs(title = paste0("H-selection with Sinkhorn’s solver"),
       x = "Number of Bins (H)", y = "OT Index") +
  theme_minimal() +
  theme(
    axis.text    = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    plot.title   = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 16, face = "bold")
  )+
  scale_color_discrete(
    labels = c(
      alpha_value   = "Covariance choice",
      form_idx      = "Covariate choice",
      mesh_idx      = "Mesh choice",
      prior_range   = "Range prior",
      prior_sigma   = "SD prior"
    )
  )+
  theme(strip.text = element_text(face = "bold"))
ggsave("Stage1_M_selection_sinkhorn_on_log_mean.png", width = 12, height = 8, dpi = 300)


message("\nStage 1 done. Inspect Stage1_M_selection.png, set H_FINAL, ",
        "then run Stage 2.")


#final run with selected binning
EPS_FINAL       <- 0.01
N_REPS_FINAL    <- 200      # 95th pct needs more than 200 to be stable
H_FINAL         <- 10       # <-- SET THIS AFTER READING THE STAGE 1 PLOT

sinkhorn_final  <- list(epsilon = EPS_FINAL,  numIterations = SINK_ITERS,
                        maxErr = SINK_MAXERR)
OUT_TAG      <- "all"

stopifnot(exists("uncertainty_factors"), exists("BETA"),
          exists("dummy_threshold"), exists("sinkhorn_final"))

# CI extractor that is ROBUST to confint() row ordering: returns a data.frame
# keyed by factor name, so we can join rather than assume positional order.
ci_by_factor <- function(obj, factor_names) {
  ci <- tryCatch(confint(obj), error = function(e) NULL)
  if (is.null(ci)) {
    # fall back: no CI available -> NA columns
    return(data.frame(Factor = factor_names,
                      Low_CI = NA_real_, High_CI = NA_real_,
                      stringsAsFactors = FALSE))
  }
  ci <- as.data.frame(ci)
  # find the id column (input names) and the low/high columns, tolerant of
  # naming differences across gsaot versions.
  id_col  <- intersect(c("input", "Input", "variable", "term"), names(ci))
  low_col <- intersect(c("low.ci", "conf.low", "lower", "low"), names(ci))
  hi_col  <- intersect(c("high.ci", "conf.high", "upper", "high"), names(ci))
  if (length(id_col) == 0) {
    # no id column: assume row order matches factor_names (last resort)
    ci$Factor <- factor_names
    id_col <- "Factor"
  }
  data.frame(
    Factor  = as.character(ci[[id_col[1]]]),
    Low_CI  = if (length(low_col)) ci[[low_col[1]]] else NA_real_,
    High_CI = if (length(hi_col))  ci[[hi_col[1]]]  else NA_real_,
    stringsAsFactors = FALSE
  )
}

run_stage2_for_M <- function(M_FINAL) {
  message("\n=== STAGE 2: M = ", M_FINAL, " ===")
  stopifnot(N / M_FINAL >= 100)
  
  # Global multivariate index (bootstrapped) + CI + dummy floor
  glob_obj       <- ot_indices(x = uncertainty_factors, y = BETA, M = M_FINAL,
                               solver = "sinkhorn", boot = TRUE, R = BOOT_R,
                               solver_optns = sinkhorn_final)
  glob_idx_final <- glob_obj$indices
  glob_ci        <- ci_by_factor(glob_obj, colnames(uncertainty_factors))
  
  set.seed(SEED)
  glob_thr_final <- dummy_threshold(BETA, M_FINAL, "sinkhorn",
                                    sinkhorn_final, N_REPS_FINAL, SEED)
  
  # Per-PC exact-1d indices (bootstrapped) + CI + dummy floor
  pc_final <- lapply(1:n_pcs_to_test, function(i) {
    obj <- ot_indices_1d(x = uncertainty_factors, y = BETA[, i],
                         M = M_FINAL, boot = TRUE, R = BOOT_R)
    thr <- dummy_threshold(BETA[, i], M_FINAL, "1d", NULL,
                           N_REPS_FINAL, SEED + i)
    list(index = obj$indices,
         ci    = ci_by_factor(obj, colnames(uncertainty_factors)),
         threshold = thr)
  })
  
  # Assemble long results, joining CIs BY FACTOR NAME within each block.
  block <- function(idx_vec, ci_df, thr, type_label) {
    d <- data.frame(Factor = colnames(uncertainty_factors),
                    Index  = as.numeric(idx_vec),
                    Threshold = thr,
                    Type = type_label,
                    stringsAsFactors = FALSE)
    left_join(d, ci_df, by = "Factor")
  }
  
  parts <- c(
    list(block(glob_idx_final, glob_ci, glob_thr_final, "Global")),
    lapply(1:n_pcs_to_test, function(i)
      block(pc_final[[i]]$index, pc_final[[i]]$ci,
            pc_final[[i]]$threshold, paste0("PC", i)))
  )
  final_results <- bind_rows(parts) %>%
    mutate(M = M_FINAL,
           Relevant        = Index > Threshold,
           Signif_vs_floor = Low_CI > Threshold)
  
  write.csv(final_results,
            sprintf("Final_results_sinkhorn_more_M%d_%s.csv", M_FINAL, OUT_TAG),
            row.names = FALSE)
  print(final_results)
  
  # Figure with CI error bars + irrelevance floor
  final_band <- final_results %>% group_by(Type) %>%
    summarise(Threshold = mean(Threshold), .groups = "drop")
  
  p <- ggplot(final_results, aes(Factor, Index, fill = Relevant)) +
    geom_col() +
    geom_errorbar(aes(ymin = Low_CI, ymax = High_CI), width = 0.3) +
    geom_hline(data = final_band, aes(yintercept = Threshold),
               color = "firebrick4", linetype = "dotted", linewidth = 1.8) +
    facet_wrap(~Type) + coord_cartesian(ylim = c(0, 1))+
    scale_fill_manual(values = c(`TRUE` = "steelblue4", `FALSE` = "grey70")) +
    labs(title = paste0("Which factors drive variation in the scenario scores? (H = ", M_FINAL,")"),
         x = NULL, y = "OT Index") +
    theme_minimal() +
    theme(
      axis.text    = element_text(size = 18, face = "bold"),
      axis.title   = element_text(size = 18, face = "bold"),
      legend.text  = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 18, face = "bold"),
      plot.title   = element_text(size = 18, face = "bold"),
      strip.text = element_text(size = 16, face = "bold")
    )+  scale_x_discrete(
      labels = c(
        alpha_value   = "Covariance choice",
        form_idx      = "Covariate choice",
        mesh_idx      = "Mesh choice",
        prior_range   = "Range prior",
        prior_sigma   = "SD prior"
      )
    )+
    theme(strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(sprintf("Final_indices_sinkhorn_more_M%d_%s.png", M_FINAL, OUT_TAG),
         p, width = 20, height = 12, dpi = 300)
  
  final_results
}

#individual plots for each binning
all_M_results <- lapply(M_LIST, run_stage2_for_M) %>% bind_rows()

write.csv(all_M_results,
          sprintf("Final_results_sinkhorn_allM_%s.csv", OUT_TAG), row.names = FALSE)

# Combined comparison across M (Global rows), to see if the verdict is M-stable
glob_all <- all_M_results %>% filter(Type == "Global")
ggplot(glob_all, aes(Factor, Index, fill = Signif_vs_floor)) +
  geom_col() +
  geom_errorbar(aes(ymin = Low_CI, ymax = High_CI), width = 0.3) +
  geom_hline(aes(yintercept = Threshold),
             color = "firebrick4", linetype = "dotted", linewidth = 0.8) +
  facet_wrap(~paste0("M=", M), scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "steelblue4", `FALSE` = "grey70"),
                    na.value = "grey85") +
  labs(title = paste0("Which factors drive variation in the scenario scores? (H = ", M_FINAL,")"),
       x = NULL, y = "OT Index")+
  theme_minimal() +
  theme(
    axis.text    = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    plot.title   = element_text(size = 18, face = "bold"),
    plot.subtitle= element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold")
  )+  scale_x_discrete(
    labels = c(
      alpha_value   = "Covariance choice",
      form_idx      = "Covariate choice",
      mesh_idx      = "Mesh choice",
      prior_range   = "Range prior",
      prior_sigma   = "SD prior"
    )
  )+
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf("Final_indices_sinkhorn_allM_%s.png", OUT_TAG),
       width = 15, height = 8, dpi = 300)



ggplot(glob_all, aes(Factor, Index, fill = Signif_vs_floor)) +
  geom_col() +
  geom_errorbar(aes(ymin = Low_CI, ymax = High_CI), width = 0.3) +
  geom_hline(aes(yintercept = Threshold),
             color = "firebrick4", linetype = "dotted", linewidth = 0.8) +
  facet_wrap(~paste0("M=", M), scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "steelblue4", `FALSE` = "grey70"),
                    na.value = "grey85") +
  labs(title = paste0("OT sensitivity indices (H = ", M_FINAL,")"),
       x = NULL, y = "OT Index")+
  theme_minimal() +
  theme(
    axis.text    = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    plot.title   = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 16, face = "bold")
  )+  scale_color_discrete(
    labels = c(
      alpha_value   = "Covariance choice",
      form_idx      = "Covariate choice",
      mesh_idx      = "Mesh choice",
      prior_range   = "Range prior",
      prior_sigma   = "SD prior",
      retain_metric = "Filtering criterion"
    )
  )+
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf("Final_indices_sinkhorn_allM_%_nosubtitles.png", OUT_TAG),
       width = 13, height = 6, dpi = 300)

message("\nWrote per-M CSVs/PNGs, Final_results_allM_", OUT_TAG,
        ".csv and the across-M comparison figure.")
message("Signif_vs_floor = Low_CI > Threshold is the defensible relevance ",
        "test (CI lower bound clears the irrelevance floor).")



#################### reload result ############
library(readr)
final_results <- read_csv("Final_results_sinkhorn_more_M10_all.csv")

final_band <- final_results %>%  group_by(Type) %>%  summarise(Threshold = mean(Threshold),    .groups = "drop" )


final<-ggplot(final_results, aes(Factor, Index, fill = Relevant)) +
  geom_col() +
  geom_errorbar(
    aes(ymin = Low_CI, ymax = High_CI),
    width = 0.3
  ) +
  geom_hline(
    data = final_band,
    aes(yintercept = Threshold),
    color = "firebrick4",
    linetype = "dotted",
    linewidth = 1.8
  ) +
  facet_wrap(~Type) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(
    values = c(`TRUE` = "steelblue4",
               `FALSE` = "grey70")
  ) +
  scale_x_discrete(
    labels = c(
      alpha_value = "Covariance choice",
      form_idx    = "Covariate choice",
      mesh_idx    = "Mesh choice",
      prior_range = "Range prior",
      prior_sigma = "SD prior"
    )
  ) +  labs(title = paste0("Which factors drive variation in the scenario scores? (H = ", M_FINAL,")"),
            x = NULL, y = "OT Index")+
  theme_minimal() +
  theme(
    axis.text    = element_text(size = 26, face = "bold"),
    axis.title   = element_text(size = 26, face = "bold"),
    legend.text  = element_text(size = 26, face = "bold"),
    legend.title = element_text(size = 26, face = "bold"),
    plot.title   = element_text(size = 30, face = "bold"),
    #plot.subtitle= element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 24, face = "bold")
  )+    theme(strip.text = element_text(face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf("Final_indices_sinkhorn_more_M%d_%s.png", 10, "all"),
       final, width = 22, height = 15, dpi = 300)

