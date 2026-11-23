###############################################################################
# Optimal-transport sensitivity analysis (OT-SA) of the retained INLA ensemble
#
#   1.  H-selection: Sinkhorn OT indices and irrelevance thresholds over a grid
#            of bin numbers H (global and per PC)
#   2.  final analysis at H_FINAL:
#              - Sinkhorn OT indices, global and per PC, with bootstrap CIs
#              - irrelevance thresholds from uniform dummy inputs
#              - exact OT indices (network simplex globally, exact 1D per PC)
#              - Wasserstein-Bures (WB) decomposition: advective + diffusive,
#                non-Gaussian = exact OT - (advective + diffusive),
#                regularisation bias = Sinkhorn - exact
#   Figures  H-selection; indices with CI and threshold; Sinkhorn index | WB decomposition
#            remainder = Sinkhorn - (advective + diffusive) = non-Gaussian + regularisation bias
#
# Author: Corinna Perchtold
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gsaot)
library(future.apply)

stopifnot(packageVersion("gsaot") >= "1.1.1")

# =============================================================================
# Settings
# =============================================================================
UNCERTAINTY_FILE <- "01_Uncertainty_18k_preprocessed.R"
PCA_FILE         <- "01_PCA_RMSE_log10_Cov95.R"
PCA_OBJ          <- "pca_results_mean_rmse_log10_cov95"
METRICS          <- c("RMSE_log10", "Coverage95")      # retention criteria (union)
OUT_TAG          <- "RMSE_log10_Cov95"

FACTORS     <- c("form_idx", "mesh_idx", "alpha_value", "prior_range", "prior_sigma")
CATEGORICAL <- c("form_idx", "mesh_idx", "alpha_value")  # priors stay numeric
FACTOR_LABELS <- c(form_idx = "Covariates", mesh_idx = "Mesh", alpha_value = "Covariance",
                   prior_range = "Range prior", prior_sigma = "SD prior")

SEED            <- 20260624
THRESH_QUANTILE <- 0.95                  # upper edge of the irrelevance zone
SINKHORN        <- list(epsilon = 0.01, numIterations = 5e5, maxErr = 1e-4)
MIN_PER_BIN     <- 100                   # minimum scenarios per bin

H_GRID        <- c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25)   # Stage 1
N_REPS_SELECT <- 50                      # dummy draws per threshold, Stage 1
H_FINAL       <- 10                      # set after inspecting the Stage 1 figure
N_REPS_FINAL  <- 200                     # dummy draws per threshold, Stage 2
BOOT_R        <- 200                     # bootstrap replicates, Stage 2

RUN_STAGE1 <- TRUE
FS <- 26                                 # figure font size

plan(multisession)                       # set workers = n to cap the number of cores


# =============================================================================
# Data: design matrix of the modelling choices and PC scores
# =============================================================================
load(UNCERTAINTY_FILE)                   # labelled
e <- new.env(); load(PCA_FILE, envir = e)
pca <- get(PCA_OBJ, envir = e); rm(e)

labelled <- labelled %>% filter(retain_metric %in% METRICS)
design   <- labelled[, FACTORS] %>% mutate(across(all_of(CATEGORICAL), as.factor))

K    <- pca$n_pcs
BETA <- as.matrix(pca$BETA)[, 1:K, drop = FALSE]
if (!is.null(rownames(BETA)) && "scenario_id" %in% names(labelled)) {   # align by row name
  key <- paste0(labelled$scenario_id, "_", labelled$retain_metric)
  if (all(key %in% rownames(BETA))) BETA <- BETA[key, , drop = FALSE]
}
stopifnot(nrow(BETA) == nrow(design))

N     <- nrow(design)
TYPES <- c("Global", paste0("PC", 1:K))
stopifnot(N / max(H_GRID) >= MIN_PER_BIN, N / H_FINAL >= MIN_PER_BIN)
message(sprintf("N = %d scenarios, K = %d PCs, %d factors; %.0f scenarios per bin at H = %d",
                N, K, length(FACTORS), N / H_FINAL, H_FINAL))


# =============================================================================
# Helpers
# =============================================================================
# irrelevance threshold: 95th percentile of the index of a uniform dummy input
dummy_threshold <- function(y, H, n_reps, seed) {
  draws <- future_replicate(n_reps, {
    irrelevance_threshold(y = y, M = H, solver = "sinkhorn", solver_optns = SINKHORN,
                          dummy_optns = list(distr = "runif"))$indices
  }, future.seed = seed)
  as.numeric(quantile(draws, THRESH_QUANTILE, names = FALSE))
}

# y for block "Global" (all PCs) or "PCj" (one PC)
y_of <- function(type) if (type == "Global") BETA else BETA[, as.integer(sub("PC", "", type)), drop = FALSE]

# bootstrap CI by factor name (robust to the row order of confint())
ci_by_factor <- function(obj) {
  ci <- tryCatch(as.data.frame(confint(obj)), error = function(e) NULL)
  if (is.null(ci)) return(data.frame(Factor = FACTORS, Low_CI = NA_real_, High_CI = NA_real_))
  id <- intersect(c("input", "Input", "variable"), names(ci))
  if (length(id) == 0) { ci$input <- FACTORS; id <- "input" }
  data.frame(Factor = as.character(ci[[id[1]]]), Low_CI = ci$low.ci, High_CI = ci$high.ci)
}

theme_ot <- theme_minimal() +
  theme(text         = element_text(size = FS, face = "bold"),
        axis.title   = element_text(size = FS, face = "bold"),
        axis.text    = element_text(size = FS, face = "bold", colour = "black"),
        strip.text   = element_text(size = FS, face = "bold"),
        legend.text  = element_text(size = FS, face = "bold"),
        legend.title = element_text(size = FS, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1.5, "lines"))


# =============================================================================
# Stage 1: H-selection (Sinkhorn, global and per PC, no bootstrap)
# =============================================================================
if (RUN_STAGE1) {
  message("\n=== Stage 1: H-selection ===")
  stage1 <- bind_rows(lapply(H_GRID, function(H) {
    message("  H = ", H)
    bind_rows(lapply(seq_along(TYPES), function(t) {
      y <- y_of(TYPES[t])
      data.frame(H = H, Type = TYPES[t], Factor = FACTORS,
                 Index = as.numeric(ot_indices(x = design, y = y, M = H, solver = "sinkhorn",
                                               solver_optns = SINKHORN)$indices),
                 Threshold = dummy_threshold(y, H, N_REPS_SELECT, SEED + t))
    }))
  }))
  write.csv(stage1, sprintf("Stage1_H_selection_%s.csv", OUT_TAG), row.names = FALSE)

  band1 <- distinct(stage1, H, Type, Threshold) %>% mutate(Type = factor(Type, levels = TYPES))
  p1 <- ggplot(mutate(stage1, Type = factor(Type, levels = TYPES))) +
    geom_ribbon(data = band1, aes(H, ymin = 0, ymax = Threshold), fill = "firebrick", alpha = 0.15) +
    geom_line(data = band1, aes(H, Threshold), colour = "firebrick4", linetype = "dotted", linewidth = 1.4) +
    geom_line(aes(H, Index, colour = Factor), linewidth = 1.6) +
    geom_point(aes(H, Index, colour = Factor), size = 4) +
    geom_vline(xintercept = H_FINAL, colour = "grey40", linetype = "dashed", linewidth = 1) +
    facet_wrap(~ Type, ncol = 2) + coord_cartesian(ylim = c(0, 1)) +
    scale_colour_discrete(labels = FACTOR_LABELS) +
    labs(x = "Number of bins (H)", y = "OT index", colour = NULL) + theme_ot +
    guides(colour = guide_legend(nrow = 2))
  ggsave(sprintf("Stage1_H_selection_%s.png", OUT_TAG), p1,
         width = 14, height = 6.5 * ceiling(length(TYPES) / 2) + 2.5, dpi = 300, bg = "white")
}


# =============================================================================
# Stage 2: final analysis at H_FINAL
# =============================================================================
message("\n=== Stage 2: H = ", H_FINAL, " ===")
H <- H_FINAL

final <- bind_rows(lapply(seq_along(TYPES), function(t) {
  type <- TYPES[t]; y <- y_of(type)
  message("  ", type)

  # Sinkhorn index with bootstrap CI (bootstrap indices are bias-corrected by gsaot)
  sink <- ot_indices(x = design, y = y, M = H, solver = "sinkhorn",
                     solver_optns = SINKHORN, boot = TRUE, R = BOOT_R)
  thr  <- dummy_threshold(y, H, N_REPS_FINAL, SEED + t)

  # exact OT: network simplex for the multivariate scores, exact 1D solver per PC
  exact <- if (type == "Global") ot_indices(x = design, y = y, M = H, solver = "transport")$indices
           else ot_indices_1d(x = design, y = as.numeric(y), M = H)$indices

  # WB decomposition (Gaussian approximation): advective + diffusive
  wb <- ot_indices_wb(x = design, y = y, M = H)

  data.frame(Type = type, Factor = FACTORS,
             Index = as.numeric(sink$indices), Threshold = thr) %>%
    left_join(ci_by_factor(sink), by = "Factor") %>%
    mutate(Exact     = as.numeric(exact),
           WB        = as.numeric(wb$indices),
           Advective = as.numeric(unlist(wb$adv)),
           Diffusive = as.numeric(unlist(wb$diff)))
})) %>%
  mutate(NonGaussian     = pmax(Exact - (Advective + Diffusive), 0),   # exact OT - WB
         Reg_bias        = Index - Exact,                              # Sinkhorn - exact
         Remainder       = Index - (Advective + Diffusive),            # grey part of the figure
         Relevant        = Index > Threshold,
         Signif_vs_floor = Low_CI > Threshold,
         H               = H)

saveRDS(final, sprintf("Final_results_sinkhornPC_M%d.rds", H))
write.csv(final, sprintf("Final_results_sinkhornPC_M%d_%s.csv", H, OUT_TAG), row.names = FALSE)
print(final %>% mutate(across(where(is.numeric), ~ round(.x, 3))))


# =============================================================================
# Figures
# =============================================================================
res    <- final %>% mutate(Type = factor(Type, levels = TYPES), Factor = factor(Factor, levels = FACTORS))
thr_df <- distinct(res, Type, Threshold)
fig_h  <- 6.5 * ceiling(length(TYPES) / 2) + 2.5

# (1) Sinkhorn indices with 95% CI and irrelevance threshold
p_idx <- ggplot(res, aes(Factor, Index, fill = Signif_vs_floor)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Low_CI, ymax = High_CI), width = 0.3, linewidth = 1.2) +
  geom_hline(data = thr_df, aes(yintercept = Threshold), colour = "firebrick4",
             linetype = "dotted", linewidth = 1.8) +
  facet_wrap(~ Type, ncol = 2) + coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c(`TRUE` = "steelblue4", `FALSE` = "grey70"),
                    labels = c(`TRUE` = "yes", `FALSE` = "no"), name = "Relevant") +
  scale_x_discrete(labels = FACTOR_LABELS) +
  labs(x = NULL, y = "OT index") + theme_ot +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())
ggsave(sprintf("OT_indices_%s_H%d.png", OUT_TAG, H), p_idx,
       width = 14, height = fig_h, dpi = 300, bg = "white")

# (2) per factor: Sinkhorn index (CI) | advective + diffusive + remainder up to the Sinkhorn index
if (any(res$Remainder < 0, na.rm = TRUE))
  message(sum(res$Remainder < 0, na.rm = TRUE),
          " factor(s) with Sinkhorn index below advective + diffusive: remainder set to 0")
res <- res %>% mutate(Remainder = pmax(Remainder, 0))
fx  <- setNames(seq_along(FACTORS), FACTORS)
off <- 0.2
bars <- bind_rows(
  res %>% transmute(Type, Factor, Component = "Sinkhorn", value = Index,
                    x = fx[as.character(Factor)] - off),
  res %>% dplyr::select(Type, Factor, Advective, Diffusive, Remainder) %>%
    pivot_longer(c(Advective, Diffusive, Remainder), names_to = "Component", values_to = "value") %>%
    mutate(x = fx[as.character(Factor)] + off)) %>%
  mutate(Component = factor(Component,
                            levels = c("Sinkhorn", "Advective", "Diffusive", "Remainder"),
                            labels = c("OT index (Sinkhorn)", "Advective", "Diffusive",
                                       "Non-Gaussian and regularisation bias")))
ci_df <- res %>% transmute(Type, x = fx[as.character(Factor)] - off, Low_CI, High_CI)

p_dec <- ggplot(bars, aes(x, value, fill = Component)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.38) +
  geom_errorbar(data = ci_df, aes(x = x, ymin = Low_CI, ymax = High_CI), inherit.aes = FALSE,
                width = 0.15, linewidth = 1.2) +
  geom_hline(data = thr_df, aes(yintercept = Threshold), colour = "firebrick4",
             linetype = "dotted", linewidth = 1.8) +
  facet_wrap(~ Type, ncol = 2) + coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = fx, labels = FACTOR_LABELS) +
  scale_fill_manual(values = c("steelblue4", "#6baed6", "#fd8d3c", "grey75")) +
  labs(x = NULL, y = "OT index", fill = NULL) + theme_ot +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank()) +
  guides(fill = guide_legend(nrow = 2))
ggsave(sprintf("OT_indices_WB_%s_H%d.png", OUT_TAG, H), p_dec,
       width = 14, height = fig_h, dpi = 300, bg = "white")

# check: the WB index is the sum of its two terms, and exact >= WB (Gelbrich bound)
chk <- final %>% summarise(max_wb_gap   = max(abs(WB - Advective - Diffusive)),
                           n_exact_lt_wb = sum(Exact < WB - 1e-8))
print(chk)
message("Done. Signif_vs_floor = lower CI above the irrelevance threshold.")
