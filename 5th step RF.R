###############################################################################
# Second spatial engine: the same workflow with a random forest (RF)
#
#   1. Uncertainty matrix for RF (all factors discrete -> partition by level, no binning)
#   2. Spatial cross-validation of every scenario (RMSE, R2 on log10; 95% coverage
#      from a quantile regression forest)
#   3. Retain the best 10% by RMSE  U  best 10% by 95% coverage (same rule as INLA)
#   4. Refit the retained scenarios on all stations and predict the grid (log10)
#   5. PCA of the map ensemble (centred, scaled, 70% variance), as for INLA
#   6. OT-SA with the same settings as the INLA analysis:
#        Sinkhorn (eps = 0.01) for the global map AND for each PC, bootstrap R = 200,
#        uniform-dummy thresholds (95th percentile of 200 draws, H = 10)
#   7. WB decomposition (advective, diffusive, non-Gaussian = exact - WB;
#      exact = transport solver for the global map, 1D quantile formula per PC) and 2x2 figure
#
# Every stage is saved; rerunning the script skips finished stages.
# Same data as the INLA prediction step (log10 concentrations, same covariates, same grid).
#
# Author : Corinna Perchtold
###############################################################################

library(ranger)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gsaot)
library(future)
library(future.apply)
library(progressr)

library(sf)
stopifnot(as.character(packageVersion("gsaot")) == "1.1.1")   



load("01_spatial_points.RData")             
load("01_Uncertainty_18k_preprocessed.RData")   

# ---- settings ----------------------------------------------------------------
SEED        <- 20261103
OUT_DIR     <- "./rf_workflow"
N_TREES     <- 1000
USE_OPTIONAL <- TRUE                   # add sample fraction and split rule as factors
N_SEEDS     <- 3                       # seed repetitions as a factor (1 = no seed factor)
N_FOLDS     <- 10                      # spatial blocks (k-means on station coordinates)
Q_RETAIN    <- 0.10                    # as for INLA
VAR_THRESH  <- 0.70                    # as for INLA
H           <- 10                      # only used by the uniform dummy (all RF factors are discrete)
BOOT_R      <- 200
N_REPS_THR  <- 200
THRESH_Q    <- 0.95
SINK_OPTS   <- list(epsilon = 0.01, numIterations = 5e5, maxErr = 1e-4)
N_WORKERS   <- max(1, min(32, parallelly::availableCores() - 1))
MIN_PER_LEVEL <- 100                   # guideline per partition set (Chiani et al.)
dir.create(OUT_DIR, showWarnings = FALSE)
f_out <- function(...) file.path(OUT_DIR, ...)

stations_sf_113AK <- as.data.frame(coords.stations2) %>%
  st_as_sf(coords = c(1, 2), crs = 2154)

# same data as the INLA prediction step
y_all  <- log10(obs.data2$resultat)
X_obs  <- covars.data2
X_grid <- spatial_points_113AK_df[, 1:6]
stopifnot(length(y_all) == nrow(coords.stations2), nrow(X_grid) == nrow(coords.prediction2))


# =============================================================================
# 1. Uncertainty matrix
# =============================================================================
# Covariate sets from the INLA formulas. Linear and RW1 terms of the same covariate
# are identical inputs for a tree model, so formulas that differ only in the functional
# form collapse to one covariate SET. The factor is therefore "covariate set".
vars_of  <- function(f) sort(setdiff(all.vars(f), c("y", "Intercept", "spatial_field", "spde")))
cov_sets_all <- lapply(formula_list, vars_of)
set_key  <- sapply(cov_sets_all, paste, collapse = "+")
cov_sets <- cov_sets_all[!duplicated(set_key)]
names(cov_sets) <- paste0("S", seq_along(cov_sets))
cat(sprintf("%d INLA formulas -> %d distinct covariate sets for RF\n", length(formula_list), length(cov_sets)))
print(data.frame(set = names(cov_sets), covariates = sapply(cov_sets, paste, collapse = ", ")))
stopifnot(all(unlist(cov_sets) %in% names(X_obs)), all(unlist(cov_sets) %in% names(X_grid)))

# Factors (all discrete -> partitioned by level):
#   cov_set   : which covariates enter                  (analogue of the INLA formula)
#   spatial   : how location enters                     (analogue of the SPDE field / mesh)
#               none   = covariates only
#               coords = + x, y coordinates
#               edf    = + x, y and distances to the 4 corners and the centre of the
#                        domain (Euclidean distance fields, Behrens et al. 2018)
#               rfsp   = + distances to every TRAINING station (Hengl et al. 2018)
#   mtry_rule : features tried per split: sqrt(p), p/3, p (bagging)
#   min_node  : minimum terminal node size (smoothness)
#   optional  : sample_fr (share of stations per tree), splitrule (variance / extratrees)
#   seed_rep  : repetition of the forest with another seed (Monte Carlo noise reference)
grid_args <- list(cov_set   = names(cov_sets),
                  spatial   = c("none", "coords", "edf", "rfsp"),
                  mtry_rule = c("sqrt", "third", "all"),
                  min_node  = c(5, 10, 15, 20))
if (USE_OPTIONAL) grid_args <- c(grid_args, list(sample_fr = c(0.5, 0.632, 1),
                                                 splitrule = c("variance", "extratrees")))
if (N_SEEDS > 1)  grid_args <- c(grid_args, list(seed_rep = seq_len(N_SEEDS)))
RF_UNC <- do.call(expand.grid, c(grid_args, stringsAsFactors = FALSE))
if (!"sample_fr" %in% names(RF_UNC)) RF_UNC$sample_fr <- 1
if (!"splitrule" %in% names(RF_UNC)) RF_UNC$splitrule <- "variance"
if (!"seed_rep"  %in% names(RF_UNC)) RF_UNC$seed_rep  <- 1
RF_UNC$scenario_id <- seq_len(nrow(RF_UNC))
fac_rf <- names(grid_args)
cat("RF scenarios:", nrow(RF_UNC), "\n")

# expected retained rows and rows per level (OT guideline: >= 100 per partition set)
n_ret <- 2 * Q_RETAIN * nrow(RF_UNC)
lev   <- sapply(RF_UNC[, fac_rf], n_distinct)
print(data.frame(levels = lev, expected_rows_per_level = round(n_ret / lev)))
if (any(n_ret / lev < MIN_PER_LEVEL))
  warning("some factors will have fewer than ", MIN_PER_LEVEL, " retained rows per level; ",
          "reduce their levels or enlarge the matrix")
saveRDS(list(RF_UNC = RF_UNC, cov_sets = cov_sets), f_out("01_rf_uncertainty_matrix.rds"))


# =============================================================================
# 2. Spatial cross-validation
# =============================================================================
# EDF reference points: corners and centre of the PREDICTION domain (projected coordinates,
# fixed, independent of the stations -> identical in every CV fold and on the grid)
ref_pts <- rbind(c(bb[1, 1], bb[1, 2]),   # (x_min, y_min)  lower-left corner
                 c(bb[1, 1], bb[2, 2]),   # (x_min, y_max)  upper-left corner
                 c(bb[2, 1], bb[1, 2]),   # (x_max, y_min)  lower-right corner
                 c(bb[2, 1], bb[2, 2]),   # (x_max, y_max)  upper-right corner
                 colMeans(bb))            # ((x_min+x_max)/2, (y_min+y_max)/2)  centre                                

dist_to <- function(xy, pts, prefix) {
  D <- sapply(seq_len(nrow(pts)), function(k) sqrt((xy[, 1] - pts[k, 1])^2 + (xy[, 2] - pts[k, 2])^2))
  D <- matrix(D, nrow = nrow(xy))
  colnames(D) <- paste0(prefix, seq_len(nrow(pts)))
  D
}

# xy     : locations of the rows (stations or grid points)
# xy_ref : TRAINING stations of the current fit (only used by "rfsp")
make_features <- function(vars, rep, X_cov, xy, xy_ref) {
  feat <- as.data.frame(X_cov[, vars, drop = FALSE])
  if (rep %in% c("coords", "edf")) { feat$X_coord <- xy[, 1]; feat$Y_coord <- xy[, 2] }
  if (rep == "edf")  feat <- cbind(feat, dist_to(xy, ref_pts, "EDF"))      # Behrens et al. (2018)
  if (rep == "rfsp") feat <- cbind(feat, dist_to(xy, xy_ref,  "dist_st"))  # Hengl et al. (2018)
  if (ncol(feat) == 0) stop("no features: empty covariate set without spatial representation")
  feat
}

mtry_of <- function(rule, p) switch(rule,
                                    sqrt  = max(1, floor(sqrt(p))),       # ranger default for regression
                                    third = max(1, round(p / 3)),         # randomForest default for regression
                                    all   = p)                            # bagging

fit_rf <- function(feat, y, row, seed) {
  ranger(x = feat, y = y, num.trees = N_TREES,
         mtry            = mtry_of(row$mtry_rule, ncol(feat)),
         min.node.size   = row$min_node,
         sample.fraction = row$sample_fr,
         splitrule       = row$splitrule,
         respect.unordered.factors = "order",
         quantreg        = TRUE, keep.inbag = FALSE,
         seed = seed, num.threads = 1)
}

# seed: differs between seed repetitions, identical across folds of one scenario
seed_of <- function(row) SEED + 100000 * row$seed_rep + row$scenario_id

# categorical covariates must have the same levels for stations and grid
for (v in intersect(names(X_obs), names(X_grid))) {
  if (is.factor(X_obs[[v]]) || is.character(X_obs[[v]]) || is.factor(X_grid[[v]])) {
    lv <- sort(unique(c(as.character(X_obs[[v]]), as.character(X_grid[[v]]))))
    X_obs[[v]]  <- factor(X_obs[[v]],  levels = lv)
    X_grid[[v]] <- factor(X_grid[[v]], levels = lv)
  }
}


predict_rf <- function(fit, feat) {
  mu <- predict(fit, data = feat, type = "response")$predictions
  q  <- predict(fit, data = feat, type = "quantiles", quantiles = c(0.025, 0.975))$predictions
  list(mean = mu, q025 = q[, 1], q975 = q[, 2])
}

set.seed(SEED)
folds <- kmeans(coords.stations2, centers = N_FOLDS, nstart = 50)$cluster
saveRDS(folds, f_out("02_rf_folds.rds"))
print(table(folds))

cv_rf_scenario <- function(i) {
  row  <- RF_UNC[i, ]
  vars <- cov_sets[[row$cov_set]]
  mu <- q025 <- q975 <- rep(NA_real_, length(y_all))
  for (k in seq_len(N_FOLDS)) {
    tr <- which(folds != k); te <- which(folds == k)
    f_tr <- make_features(vars, row$spatial, X_obs[tr, , drop = FALSE],
                          coords.stations2[tr, , drop = FALSE], coords.stations2[tr, , drop = FALSE])
    f_te <- make_features(vars, row$spatial, X_obs[te, , drop = FALSE],
                          coords.stations2[te, , drop = FALSE], coords.stations2[tr, , drop = FALSE])
    p <- predict_rf(fit_rf(f_tr, y_all[tr], row, seed = seed_of(row)), f_te)
    mu[te] <- p$mean; q025[te] <- p$q025; q975[te] <- p$q975
  }
  data.frame(scenario_id = i,
             rmse_logscale = sqrt(mean((y_all - mu)^2)),
             R2_log10      = 1 - sum((y_all - mu)^2) / sum((y_all - mean(y_all))^2),
             Coverage_95   = mean(y_all >= q025 & y_all <= q975))  # log10 is monotone: same coverage on the original scale
}

cv_file <- f_out("02_rf_cv_stats.rds")
if (!file.exists(cv_file)) {
  print(cv_rf_scenario(1))                                     # quick test before the full run
  plan(multisession, workers = N_WORKERS)
  handlers(handler_txtprogressbar(style = 3))
  cv_stats <- with_progress({
    p <- progressor(steps = nrow(RF_UNC))
    future_lapply(RF_UNC$scenario_id, function(i) {
      out <- tryCatch(cv_rf_scenario(i), error = function(e) {
        write(sprintf("RF CV %d failed: %s", i, e$message), file = f_out("rf_error_log.txt"), append = TRUE)
        NULL })
      p(); out
    }, future.seed = TRUE, future.packages = "ranger")
  })
  plan(sequential)
  saveRDS(bind_rows(cv_stats), cv_file)
}
RF_UNC <- RF_UNC %>% left_join(readRDS(cv_file), by = "scenario_id")
cat(sprintf("Best RF CV RMSE (log10): %.3f | best INLA LGO-CV RMSE: %.3f\n",
            min(RF_UNC$rmse_logscale, na.rm = TRUE), min(UNCERTAINTY$rmse_logscale, na.rm = TRUE)))
write.csv(RF_UNC, f_out("02_rf_cv_stats.csv"), row.names = FALSE)


# =============================================================================
# 3. Retain the best 10% (same rule as for INLA)
# =============================================================================
top_rmse <- RF_UNC %>% filter(rmse_logscale < quantile(rmse_logscale, Q_RETAIN, na.rm = TRUE))
top_cov  <- RF_UNC %>% filter(Coverage_95 > quantile(Coverage_95, 1 - Q_RETAIN, na.rm = TRUE))
labelled_rf <- bind_rows(mutate(top_rmse, retain_metric = "RMSE_log10"),
                         mutate(top_cov,  retain_metric = "Coverage95")) %>%
  mutate(key = paste0(scenario_id, "_", retain_metric))
cat(sprintf("Retained rows: %d (RMSE %d, coverage %d), unique scenarios: %d\n",
            nrow(labelled_rf), nrow(top_rmse), nrow(top_cov), n_distinct(labelled_rf$scenario_id)))
cat("Coverage cut-off:", quantile(RF_UNC$Coverage_95, 1 - Q_RETAIN, na.rm = TRUE),
    "(ties at a discrete coverage value can enlarge the coverage set)\n")
print(sapply(labelled_rf[, fac_rf], function(v) min(table(v))))   # smallest level per factor
saveRDS(labelled_rf, f_out("03_rf_labelled.rds"))


# =============================================================================
# 4. Refit on all stations and predict the grid
# =============================================================================
pred_file <- function(i) f_out(sprintf("rf_scenario_%05d.rds", i))
predict_rf_scenario <- function(i) {
  if (file.exists(pred_file(i))) return(invisible(NULL))
  row    <- RF_UNC[RF_UNC$scenario_id == i, ]
  vars   <- cov_sets[[row$cov_set]]
  f_all  <- make_features(vars, row$spatial, X_obs,  coords.stations2,   coords.stations2)
  f_grid <- make_features(vars, row$spatial, X_grid, coords.prediction2, coords.stations2)
  p      <- predict_rf(fit_rf(f_all, y_all, row, seed = seed_of(row)), f_grid)
  saveRDS(list(pred_mean_logscale = p$mean, pred_q025_logscale = p$q025,
               pred_q975_logscale = p$q975), pred_file(i))
  invisible(NULL)
}
ids <- unique(labelled_rf$scenario_id)
todo <- ids[!file.exists(sapply(ids, pred_file))]
cat(length(todo), "of", length(ids), "RF maps to predict\n")
if (length(todo)) {
  plan(multisession, workers = N_WORKERS)
  handlers(handler_txtprogressbar(style = 3))
  with_progress({
    p <- progressor(steps = length(todo))
    future_lapply(todo, function(i) {
      tryCatch(predict_rf_scenario(i), error = function(e)
        write(sprintf("RF prediction %d failed: %s", i, e$message), file = f_out("rf_error_log.txt"), append = TRUE))
      p()
    }, future.seed = TRUE, future.packages = "ranger")
  })
  plan(sequential)
}
stopifnot(all(file.exists(sapply(ids, pred_file))))


# =============================================================================
# 5. PCA (as for INLA: centred, scaled, 70% variance)
# =============================================================================
pca_file <- f_out("05_rf_pca.rds")
if (!file.exists(pca_file)) {
  MAP_RF <- do.call(cbind, lapply(labelled_rf$scenario_id, function(i) readRDS(pred_file(i))$pred_mean_logscale))
  colnames(MAP_RF) <- labelled_rf$key                             # pixels x (scenario, metric)
  keep_px <- apply(MAP_RF, 1, var) > 0
  cat("Constant pixels dropped:", sum(!keep_px), "\n")
  p   <- prcomp(t(MAP_RF[keep_px, ]), center = TRUE, scale. = TRUE)
  cv  <- cumsum(p$sdev^2) / sum(p$sdev^2)
  K   <- which(cv >= VAR_THRESH)[1]
  saveRDS(list(pca_results = p, n_pcs = K, BETA = p$x[, 1:K, drop = FALSE],
               BASIS = t(p$rotation[, 1:K, drop = FALSE]), keep_px = keep_px,
               variance_explained = p$sdev^2 / sum(p$sdev^2),
               MEAN = rowMeans(MAP_RF), SD = apply(MAP_RF, 1, sd)), pca_file)
}
pca_rf <- readRDS(pca_file)
K      <- pca_rf$n_pcs
cat(sprintf("RF ensemble: %d PCs explain %.1f%% (PC1-3: %s)\n", K,
            100 * sum(pca_rf$variance_explained[1:K]),
            paste(round(100 * pca_rf$variance_explained[1:3], 1), collapse = ", ")))

# design and scores, aligned by key
BETA_RF <- pca_rf$BETA[labelled_rf$key, , drop = FALSE]
design_rf <- labelled_rf[, fac_rf] %>% mutate(across(everything(), factor))
stopifnot(identical(rownames(BETA_RF), labelled_rf$key))


# =============================================================================
# 6. OT-SA (same settings as for INLA)
# =============================================================================
one_dummy <- function(y, solver, solver_optns) {
  args <- list(y = y, M = H, solver = solver, solver_optns = solver_optns, dummy_optns = list(distr = "runif"))
  if ("R_irr" %in% names(formals(gsaot::irrelevance_threshold))) args$R_irr <- 1
  out <- do.call(gsaot::irrelevance_threshold, args)
  as.numeric(if (is.list(out)) out$indices else out)
}
dummy_threshold <- function(y, solver, solver_optns, seed) {
  draws <- future_replicate(N_REPS_THR, one_dummy(y, solver, solver_optns), future.seed = seed)
  as.numeric(quantile(draws, THRESH_Q, na.rm = TRUE))
}
ci_by_factor <- function(obj, factor_names) {
  ci <- tryCatch(as.data.frame(confint(obj)), error = function(e) NULL)
  if (is.null(ci)) return(data.frame(Factor = factor_names, Low_CI = NA_real_, High_CI = NA_real_, original = NA_real_))
  if ("component" %in% names(ci)) ci <- ci[ci$component == "wass-bures", ]
  data.frame(Factor = as.character(ci$input), Low_CI = ci$low.ci, High_CI = ci$high.ci, original = ci$original)
}

ot_file <- f_out("06_rf_OT_indices_sinkhornPC.csv")       # new name: old file used the 1D solver per PC
if (!file.exists(ot_file)) {
  GSAOT_LIB <- dirname(find.package("gsaot"))
  plan(multisession, workers = N_WORKERS, rscript_libs = c(GSAOT_LIB, .libPaths()))
  stopifnot(value(future(as.character(packageVersion("gsaot")))) == "1.1.1")
  
  set.seed(SEED)
  g <- ot_indices(design_rf, BETA_RF, M = H, solver = "sinkhorn", solver_optns = SINK_OPTS,
                  boot = TRUE, R = BOOT_R)
  glob <- data.frame(Type = "Global", Factor = names(design_rf), Index = as.numeric(g$indices),
                     Threshold = dummy_threshold(BETA_RF, "sinkhorn", SINK_OPTS, SEED)) %>%
    left_join(ci_by_factor(g, names(design_rf)), by = "Factor")
  pcs <- bind_rows(lapply(seq_len(K), function(i) {
    y_i <- matrix(BETA_RF[, i], ncol = 1)                       # one PC as a 1-column matrix
    set.seed(SEED + i)
    o <- ot_indices(design_rf, y_i, M = H, solver = "sinkhorn", solver_optns = SINK_OPTS,
                    boot = TRUE, R = BOOT_R)
    data.frame(Type = paste0("PC", i), Factor = names(design_rf), Index = as.numeric(o$indices),
               Threshold = dummy_threshold(y_i, "sinkhorn", SINK_OPTS, SEED + i)) %>%
      left_join(ci_by_factor(o, names(design_rf)), by = "Factor")
  }))
  plan(sequential)
  write.csv(bind_rows(glob, pcs) %>% mutate(Relevant = Index > Threshold, Signif_vs_floor = Low_CI > Threshold),
            ot_file, row.names = FALSE)
}
ot_rf <- read.csv(ot_file)
print(ot_rf %>% mutate(across(where(is.numeric), ~ round(.x, 3))))


# =============================================================================
# 7. WB decomposition and 2x2 figure
# =============================================================================
wb_file <- f_out("07_rf_WB_exact.rds")
if (!file.exists(wb_file)) {
  wb_g <- ot_indices_wb(design_rf, BETA_RF, M = H)
  ex_g <- ot_indices(design_rf, BETA_RF, M = H, solver = "transport")
  dec  <- bind_rows(
    data.frame(Type = "Global", Factor = names(design_rf), Exact = as.numeric(ex_g$indices),
               WB = as.numeric(wb_g$indices), Advective = as.numeric(wb_g$adv), Diffusive = as.numeric(wb_g$diff)),
    bind_rows(lapply(seq_len(K), function(i) {
      w <- ot_indices_wb(design_rf, matrix(BETA_RF[, i], ncol = 1), M = H)
      data.frame(Type = paste0("PC", i), Factor = names(design_rf),
                 Exact = as.numeric(ot_indices_1d(design_rf, BETA_RF[, i], M = H)$indices),
                 WB = as.numeric(w$indices), Advective = as.numeric(w$adv), Diffusive = as.numeric(w$diff))
    })))
  saveRDS(dec, wb_file)
}
dec_rf <- ot_rf %>% left_join(readRDS(wb_file), by = c("Type", "Factor")) %>%
  mutate(NonGaussian   = pmax(Exact - WB, 0),
         bias_sinkhorn = original - Exact)                      # entropic bias (plain Sinkhorn estimates)
write.csv(dec_rf, f_out("07_rf_decomposition.csv"), row.names = FALSE)

factor_labels_rf <- c(cov_set = "Covariate set", spatial = "Spatial representation",
                      mtry_rule = "mtry", min_node = "Min. node size",
                      sample_fr = "Sample fraction", splitrule = "Split rule",
                      seed_rep = "Seed")[fac_rf]
types <- c("Global", paste0("PC", seq_len(K)))
fx    <- setNames(seq_along(factor_labels_rf), names(factor_labels_rf))
off   <- 0.2
bars <- bind_rows(
  dec_rf %>% transmute(Type, Factor, Component = "Index", value = Index, x = fx[Factor] - off),
  dec_rf %>% dplyr::select(Type, Factor, Advective, Diffusive, NonGaussian) %>%
    pivot_longer(c(Advective, Diffusive, NonGaussian), names_to = "Component", values_to = "value") %>%
    mutate(x = fx[Factor] + off)) %>%
  mutate(Type = factor(Type, levels = types),
         Component = factor(Component, levels = c("Index", "Advective", "Diffusive", "NonGaussian"),
                            labels = c("OT index (Sinkhorn)", "Advective (mean shift)",
                                       "Diffusive (covariance change)", "Non-Gaussian")))
ci_df  <- dec_rf %>% transmute(Type = factor(Type, levels = types), x = fx[Factor] - off, Low_CI, High_CI)
thr_df <- dec_rf %>% distinct(Type, Threshold) %>% mutate(Type = factor(Type, levels = types))

p_rf <- ggplot(bars, aes(x, value, fill = Component)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.38) +
  geom_errorbar(data = ci_df, aes(x = x, ymin = Low_CI, ymax = High_CI), inherit.aes = FALSE, width = 0.15, linewidth = 0.8) +
  geom_hline(data = thr_df, aes(yintercept = Threshold), colour = "firebrick4", linetype = "dotted", linewidth = 1.2) +
  facet_wrap(~ Type, ncol = 2) +
  scale_x_continuous(breaks = fx, labels = factor_labels_rf) +
  scale_fill_manual(values = c("steelblue4", "#6baed6", "#fd8d3c", "grey75")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "OT index", fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title  = element_text(size = 16, face = "bold"),
        strip.text  = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        panel.grid.major.x = element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))
ggsave(f_out("RF_OT_decomposition_2x2_sinkhornPC.png"), p_rf, width = 20, height = 13, dpi = 300, bg = "white")

cat("Done.\n")


########### get PC maps
library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(sf)

OUT_DIR    <- "./rf_workflow"
L93_OFFSET <- c(703064, 6851896)                              # centred -> Lambert 93
pca_rf <- readRDS(file.path(OUT_DIR, "05_rf_pca.rds"))
K      <- pca_rf$n_pcs
keep   <- pca_rf$keep_px

# ---- coordinates in Lambert 93 ----------------------------------------------
grid_xy <- as.matrix(coords.prediction2[, 1:2])
if (max(abs(grid_xy[, 2])) < 1e6) grid_xy <- sweep(grid_xy, 2, L93_OFFSET, "+")   # centred -> L93
st_xy <- as.matrix(coords.stations2[, 1:2])
if (max(abs(st_xy[, 2])) < 1e6) st_xy <- sweep(st_xy, 2, L93_OFFSET, "+")
xy <- data.frame(x = grid_xy[, 1], y = grid_xy[, 2])
cell <- min(diff(sort(unique(xy$x))))                         # grid cell size (m)

stations_sf <- st_as_sf(data.frame(x = st_xy[, 1], y = st_xy[, 2]), coords = c("x", "y"), crs = 2154)

# aquifer outline from the mesh interior boundary (Lambert 93)
seg <- mesh_list[[1]]$segm$int
pts <- if (!is.null(seg$loc) && nrow(seg$loc) >= max(seg$idx)) seg$loc else mesh_list[[1]]$loc
aquifer_sf <- st_sfc(lapply(seq_len(nrow(seg$idx)), function(i) st_linestring(pts[seg$idx[i, ], 1:2])),
                     crs = 2154)

# ---- common map style --------------------------------------------------------
theme_map <- theme_minimal() +
  theme(plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 20, face = "bold"),
        axis.text    = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text  = element_text(size = 18, face = "bold"),
        legend.key.height = unit(1.5, "cm"))

map_layers <- function(show_stations = TRUE) list(
  geom_sf(data = aquifer_sf, colour = "black", linewidth = 0.6, inherit.aes = FALSE),
  if (show_stations) geom_sf(data = stations_sf, colour = "grey25", size = 1.8, inherit.aes = FALSE),
  coord_sf(crs = 2154, default_crs = 2154, datum = 4326, expand = FALSE),
  theme_map)

library(dplyr); library(ggplot2); library(patchwork); library(sf)

heat <- c("#ffff00", "#ffa500", "#ff4500", "#e31a1c", "#bd0026")

base_map <- function(df, fill_var) list(
  geom_raster(data = df, aes(lon, lat, fill = .data[[fill_var]])),
  geom_sf(data = aquifer_ll,  colour = "black", linewidth = 0.6, inherit.aes = FALSE),
  geom_sf(data = stations_ll, colour = "grey25", size = 1.5,    inherit.aes = FALSE),
  scale_x_continuous(
    breaks = c(1.3, 1.8, 2.3, 2.8, 3.3, 3.8),
    labels = c("1.3", "1.8", "2.3", "2.8", "3.3", "3.8")
  ),
  scale_y_continuous(
    breaks = c(48.1, 48.6, 49.1, 49.6),
    labels = c("48.1", "48.6", "49.1", "49.6")
  ),
  coord_sf(crs = 4326, xlim = c(1.3, 4.3), ylim = c(48.1, 49.65), expand = FALSE),
  labs(x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)"),
  theme_map)

# ---- mean ----------------------------------------------------------------
p_mean <- ggplot() + base_map(mean_ll, "mean_log10") +
  scale_fill_gradientn(colours = heat, name = "MEAN",
                       breaks = scales::breaks_pretty(4), labels = function(v) sprintf("%.2f", v)) +
  labs(title = "Ensemble mean on log10-scale")

# ---- CoV: more legend values ----------------------------------------------
p_cov <- ggplot() + base_map(cov_ll, "CoV") +
  scale_fill_gradientn(colours = heat, name = "CoV (%)",
                       breaks = scales::breaks_extended(n = 8), labels = function(v) sprintf("%.1f", v)) +
  labs(title = "Coefficient of variation")

# ---- PCs: variance share in the strip, legend values not rounded ------------
pc_lab <- setNames(sprintf("PC%d (%.1f%% variance)", seq_len(K), 100 * pca_rf$variance_explained[seq_len(K)]),
                   paste0("PC", seq_len(K)))
pc_ll  <- pc_ll %>% mutate(PC_lab = factor(pc_lab[sub(" .*", "", PC)], levels = pc_lab))
lim    <- max(abs(pc_ll$value), na.rm = TRUE)

p_pc <- ggplot() + base_map(pc_ll, "value") +
  facet_wrap(~ PC_lab, nrow = 1) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       limits = c(-lim, lim), na.value = "grey85", name = "Loading",
                       breaks = scales::breaks_extended(n = 7), labels = function(v) format(v, digits = 3)) +
  theme(strip.text = element_text(size = 18, face = "bold"))

# ---- one figure ------------------------------------------------------------
p_all <- (p_mean | p_cov) / p_pc + plot_layout(heights = c(1, 1))
ggsave(file.path(OUT_DIR, "RF_mean_CoV_PCs.png"), p_all,
       width = max(20, 8 * K), height = 16, dpi = 300, bg = "white")