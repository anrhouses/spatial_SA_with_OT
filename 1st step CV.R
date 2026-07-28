###############################################################################
# Uncertainty matrix definition, Cross-validation and performance filtering

#1. Load uncertainty design: model formula, mesh choice, smoothness parameter, PC priors
#2. Fit SPDE-INLA model.
#3. Perform Leave-Group-Out CV
#4. Compute statistical measures for filtering out bad performing scenarios
#5. Retain top-performing scenarios.

# Author: Corinna Perchtold
###############################################################################
library(progressr)
library(PopED)
library(ggplot2)
library(sf)
library(raster)
library(sp)
library(terra)
library(fmesher)
library(Matrix)
library(gstat)
library(future)
library(future.apply)
library(INLA)
library(dplyr)


load("01_spatial_points.R")



###################### define scenarios ########################
set.seed(1234)
load("01_uncertainty_matrix.R")



n_scenarios<-nrow(UNCERTAINTY)
cat("We have",n_scenarios, "scenarios")

#########################  Parallelisation setup ####################  
# Detect cores. Each INLA run uses PARDISO internally, so we avoid
# oversubscription by limiting parallel workers.
# Rule of thumb: total_cores / inla_threads_per_run.
# PARDISO typically uses 2-4 threads, so we set n_workers accordingly.
total_cores <- parallel::detectCores(logical = FALSE)
inla_threads <- as.integer(Sys.getenv("INLA_NUM_THREADS", unset = "2"))
n_workers <- floor(total_cores / inla_threads) 
cat(sprintf("Parallelisation: %d physical cores, %d INLA threads -> %d parallel workers\n",
            total_cores, inla_threads, n_workers))


# Force each INLA fit to respect the thread budget (prevents 64*128 explosion)
INLA::inla.setOption(num.threads = paste0(inla_threads, ":1"))

# Cap future's per-worker object transfer so it errors loudly instead of silently
options(future.globals.maxSize = 8 * 1024^3)   # 8 GB ceiling per worker

plan(multisession, workers = n_workers)

# Create results directory
results_dir <- "cv_results_checkpoints_test"
dir.create(results_dir, showWarnings = FALSE)

cat(sprintf("Starting LGOCV for %d scenarios | 1 INLA fit + inla.group.cv() each\n\n",
            n_scenarios))


######################## precompute A matrix ##########

cat("Pre-computing A matrices for all meshes...\n")
unique_meshes <- unique(UNCERTAINTY$mesh_idx)
A_matrices_by_mesh <- vector("list", length(unique_meshes))
names(A_matrices_by_mesh) <- paste0("mesh_", unique_meshes)

for (m in seq_along(unique_meshes)) {
  mesh_idx <- unique_meshes[m]
  curr_mesh <- mesh_list[[mesh_idx]]
  A_matrices_by_mesh[[m]] <- inla.spde.make.A(curr_mesh, loc = coords.stations2)
  cat(sprintf("  Pre-computed A matrix for mesh %d\n", mesh_idx))
}
cat("A matrix pre-computation complete!\n\n")



##########################   CV #########################  
# inla.group.cv() is a POST-PROCESSING function called on a fitted model.
# It automatically constructs CV groups based on the posterior correlations
# so for an SPDE model it knows exactly which
# stations are correlated and removes them together.
#
# num.level.sets controls the neighbourhood depth:
#   -1 = standard LOO-CV (no grouping, fast but ignores correlation)
#    1 = remove the single most informative neighbour per observation
#    3,5,10 = increasingly larger groups (more honest for spatial models)
#
# Output fields:
#   $mean  CV predictive mean of eta_i (here: log-scale)
#   $sd    CV predictive SD
#   $cv    CV predictive density p(y_i | y_{-g(i)})
#   $groups the automatically constructed groups
#num_level_sets <- 3  # adjust: larger = more spatial buffer


num_level_sets_grid <- c(-1, 3, 5, 10)   # -1 = LOO baseline

run_scenario_lgocv_multi <- function(i,  precomputed_A, 
                                     nls_grid = num_level_sets_grid) {
  row       <- UNCERTAINTY[i, ]
  curr_form <- formula_list[[row$form_idx]]
  mesh      <- mesh_list[[row$mesh_idx]]
  
  
  #spde      <- precomputed_spde[[i]]
  # Build ONLY the spde this scenario needs
  spde <- inla.spde2.pcmatern(
    mesh        = mesh,
    alpha       = row$alpha_value,
    prior.range = c(row$prior_range, 0.5),
    prior.sigma = c(row$prior_sigma, 0.5)
  )
  
  mesh_key  <- paste0("mesh_", row$mesh_idx)
  A_full    <- precomputed_A[[mesh_key]]
  
  # ---- Build stack ----
  stk <- inla.stack(
    data    = list(y = log_y),
    A       = list(A_full, 1),
    effects = list(inla.spde.make.index("spatial_field", spde$n.spde),
                   list(Intercept = 1, covars.data2)),
    tag = "obs")
  stk_data <- inla.stack.data(stk, spde = spde)
  obs_rows <- inla.stack.index(stk, "obs")$data
  
  # ---- Single INLA fit (reused for all CV depths) ----
  res <- inla(curr_form, family = "gaussian",
              data       = stk_data,
              quantiles  = c(0.025, 0.25, 0.5, 0.75, 0.975),
              control.compute   = list(dic = TRUE, cpo = TRUE, waic = TRUE,
                                       return.marginals.predictor = FALSE,
                                       openmp.strategy = "pardiso.parallel"),
              control.predictor = list(A = inla.stack.A(stk), compute = TRUE, link = 1),
              control.inla      = list(int.strategy = "eb", strategy = "adaptive"))
  
  # ---- Likelihood (nugget) variance: needed for correct back-transform ----
  # E[Y] = 10^(eta + 0.5 * ln(10) * sigma2_total), with
  #   sigma2_total = Var(eta)  [lgocv$sd^2]  +  sigma2_eps  [observation noise].
  # Omitting sigma2_eps under-corrects the high end -> the underestimation you saw.
  prec_marg  <- res$summary.hyperpar["Precision for the Gaussian observations", "mean"]
  sigma2_eps <- 1 / prec_marg
  
  # ---- CPO-based LPML (depth-independent; from the single fit) ----
  cpo_vals <- res$cpo$cpo[obs_rows]
  cpo_vals[cpo_vals <= 0 | is.na(cpo_vals)] <- .Machine$double.xmin
  lpml <- sum(log(cpo_vals))
  
  obs_orig <- obs.data2$resultat
  
  # ---- Loop over CV depths, reusing the single fit ----
  per_depth <- lapply(nls_grid, function(nls) {
    lgocv <- inla.group.cv(res, num.level.sets = nls)
    
    cv_log_mean <- lgocv$mean[obs_rows]
    cv_log_sd   <- lgocv$sd[obs_rows]
    cv_density  <- lgocv$cv[obs_rows]
    
    # total predictive variance for a NEW observation
    var_total <- cv_log_sd^2 + sigma2_eps
    sd_total  <- sqrt(var_total)
    
    # back-transforms
    pred_mean   <- 10^(cv_log_mean + var_total * log(10) / 2)  # arithmetic mean (corrected)
    pred_median <- 10^(cv_log_mean)                            # median (geometric mean)
    
    # observation-level quantiles (include nugget via sd_total)
    q025 <- 10^(cv_log_mean + qnorm(0.025) * sd_total)
    q250 <- 10^(cv_log_mean + qnorm(0.25)  * sd_total)
    q500 <- pred_median
    q750 <- 10^(cv_log_mean + qnorm(0.75)  * sd_total)
    q975 <- 10^(cv_log_mean + qnorm(0.975) * sd_total)
    
    # ---- Log-score from CV densities ----
    cv_density[cv_density <= 0 | is.na(cv_density)] <- .Machine$double.xmin
    log_score <- -mean(log(cv_density))
    
    # ---- Assemble per-observation results (ORIGINAL structure preserved) ----
    cv_results <- data.frame(
      obs           = obs_orig,
      log_obs       = log_y,
      mean_pred     = pred_mean,
      log_mean_pred = cv_log_mean,
      log_sd_pred   = cv_log_sd,
      sd_total      = sd_total,      # added: total predictive sd (incl. nugget)
      q025 = q025, q250 = q250, q500 = q500, q750 = q750, q975 = q975
    )
    
    # ---- Metrics on log10 scale ----
    rss_log    <- sum((cv_results$log_obs - cv_results$log_mean_pred)^2, na.rm = TRUE)
    tss_log    <- sum((cv_results$log_obs - mean(cv_results$log_obs, na.rm = TRUE))^2, na.rm = TRUE)
    r2_log10   <- 1 - (rss_log / tss_log)
    rmse_log10 <- sqrt(mean((cv_results$log_obs - cv_results$log_mean_pred)^2, na.rm = TRUE))
    
    # ---- Metrics on original scale (with log-normal mean correction) ----
    rss    <- sum((cv_results$obs - cv_results$mean_pred)^2, na.rm = TRUE)
    tss    <- sum((cv_results$obs - mean(cv_results$obs, na.rm = TRUE))^2, na.rm = TRUE)
    r2_val <- 1 - (rss / tss)
    
    rmse_orig <- sqrt(mean((cv_results$obs - cv_results$mean_pred)^2, na.rm = TRUE))
    mape      <- mean(abs((cv_results$obs - cv_results$mean_pred) / cv_results$obs), na.rm = TRUE) * 100
    
    # ---- log10-scale quantiles ----
    log_q025 <- cv_log_mean + qnorm(0.025) * sd_total
    log_q250 <- cv_log_mean + qnorm(0.25)  * sd_total
    log_q500 <- cv_log_mean                                   # median on log scale
    log_q750 <- cv_log_mean + qnorm(0.75)  * sd_total
    log_q975 <- cv_log_mean + qnorm(0.975) * sd_total

    
    
    # ---- WIS & Coverage ----
    wis_50 <- mean(scoringutils::wis(observed = cv_results$obs,
                                     predicted = as.matrix(cv_results[, c("q250","q500","q750")]),
                                     quantile_level = c(0.25, 0.5, 0.75)))
    wis_95 <- mean(scoringutils::wis(observed = cv_results$obs,
                                     predicted = as.matrix(cv_results[, c("q025","q500","q975")]),
                                     quantile_level = c(0.025, 0.5, 0.975)))
    
    coverage_95 <- mean(cv_results$obs >= cv_results$q025 &
                          cv_results$obs <= cv_results$q975, na.rm = TRUE)
    coverage_50 <- mean(cv_results$obs >= cv_results$q250 &
                          cv_results$obs <= cv_results$q750, na.rm = TRUE)
    
    stats <- data.frame(
      scenario_id    = i,
      num_level_sets = nls,
      RMSE           = rmse_orig,
      RMSE_log10     = rmse_log10,
      MAPE           = mape,
      R2             = r2_val,
      R2_log10       = r2_log10,
      WIS_50         = wis_50,
      WIS_95         = wis_95,
      Coverage_95    = coverage_95,
      Coverage_50    = coverage_50,
      LPML           = lpml,
      LogScore_LGOCV = log_score
    )
    
    # tag per-observation predictions with depth so they can be stacked
    cv_results$num_level_sets <- nls
    cv_results$scenario_id    <- i
    
    list(stats = stats, plot_data = cv_results)
  })
  
  # ---- Combine across depths ----
  stats_multi <- do.call(rbind, lapply(per_depth, `[[`, "stats"))
  preds_multi <- do.call(rbind, lapply(per_depth, `[[`, "plot_data"))
  
  # fit_stats are depth-independent (single fit)
  fit_stats <- data.frame(DIC = res$dic$dic, WAIC = res$waic$waic, LPML = lpml)
  
  rm(res, stk, stk_data); gc()
  
  return(list(
    stats     = stats_multi,   # one row per (scenario, depth)
    plot_data = preds_multi,   # per-observation, tagged by depth
    fit_stats = fit_stats
  ))
}

##########################  parallel loop with checkpoints  #########################  
# Identify which scenarios still need to be run (skip completed checkpoints)
todo_idx <- integer(0)
cv_all_stats <- vector("list", n_scenarios)
cv_all_preds <- vector("list", n_scenarios)
cv_all_fits  <- vector("list", n_scenarios)

for (i in 1:n_scenarios) {
  checkpoint_file <- file.path(results_dir, sprintf("scenario_%04d.rds", i))
  if (file.exists(checkpoint_file)) {
    result <- readRDS(checkpoint_file)
    cv_all_stats[[i]] <- result$stats
    cv_all_preds[[i]] <- result$plot_data
    cv_all_fits[[i]]  <- result$fit_stats
  } else {
    todo_idx <- c(todo_idx, i)
  }
}


if (length(todo_idx) > 0) {
  tic("parallel CV")
  
  # Run remaining scenarios in parallel with a progress bar.
  # progressr + future_lapply: each worker signals p() after finishing a scenario,
  # and the progress bar updates in the main process.
  results_list <- with_progress({

    future_lapply(todo_idx, function(i) {
      library(INLA)
      
      result <- tryCatch({
        res <- run_scenario_lgocv_multi(i,  A_matrices_by_mesh)  #spde_list,
        saveRDS(res, file.path(results_dir, sprintf("scenario_%04d.rds", i)))
        res
      }, error = function(e) {
        write(sprintf("Scenario %d failed: %s", i, e$message),
              file = "cv_error_log.txt", append = TRUE)
        NULL
      })
      
      result
    }, future.seed = TRUE,
    future.globals = c("run_scenario_lgocv_multi", "UNCERTAINTY", "formula_list",
                       "mesh_list", "A_matrices_by_mesh", 
                       "obs.data2", "covars.data2", 
                       "num_level_sets_grid", "results_dir"),
    future.packages = c("INLA", "scoringutils"))
  })
  
  # Collect results back into the master lists
  for (j in seq_along(todo_idx)) {
    i <- todo_idx[j]
    result <- results_list[[j]]
    if (!is.null(result)) {
      cv_all_stats[[i]] <- result$stats
      cv_all_preds[[i]] <- result$plot_data
      cv_all_fits[[i]]  <- result$fit_stats
    }
  }
  
  elapsed <- toc(quiet = TRUE)
  cat(sprintf("\nParallel CV done. %d scenarios in %.1f minutes (%.1fs/scenario avg)\n",
              length(todo_idx),
              (elapsed$toc - elapsed$tic) / 60,
              (elapsed$toc - elapsed$tic) / length(todo_idx)))
}

# Save final progress
saveRDS(list(
  completed = n_scenarios,
  total = n_scenarios,
  timestamp = Sys.time()
), file.path(results_dir, "progress.rds"))

# Shut down workers
plan(sequential)
cat("All done.\n")



##########################  CV filtering #########################  
library(ggplot2)
library(dplyr)

results_dir <- "~/cv_results_checkpoints"
n_scenarios <- nrow(UNCERTAINTY)

cv_all_stats <- vector("list", n_scenarios)
cv_all_preds <- vector("list", n_scenarios)
for (i in 1:n_scenarios) {
  f <- file.path(results_dir, sprintf("scenario_%04d.rds", i))
  if (file.exists(f)) {
    r <- readRDS(f)
    cv_all_stats[[i]] <- r$stats        # one row per depth
    cv_all_preds[[i]] <- r$plot_data
  }
}

# master table: rows = scenario x num_level_sets, all metrics as columns
all_stats_df <- bind_rows(cv_all_stats)

plot_metric_by_depth <- function(df, metric, bins = 80) {
  ggplot(df, aes(x = .data[[metric]], fill = factor(num_level_sets))) +
    geom_histogram(bins = bins, alpha = 0.5, position = "identity") +
    labs(title = sprintf("%s across CV depths", metric),
         x = metric, y = "count", fill = "num_level_sets") +
    theme_minimal(base_size = 14)
}
plot_metric_by_depth(all_stats_df, "RMSE_log10")
plot_metric_by_depth(all_stats_df, metric = "R2_log10")
plot_metric_by_depth(all_stats_df, metric = "MAPE")

plot_metric_by_depth(all_stats_df, metric = "Coverage_50")
plot_metric_by_depth(all_stats_df, metric = "Coverage_95")

plot_metric_by_depth(all_stats_df, metric = "WIS_95")
plot_metric_by_depth(all_stats_df, metric = "WIS_50")


####################### plot by num_level_sets #########################  
# ---- histogram of one metric at one CV depth ----
# change `metric` and `depth` to whatever you want to inspect
plot_metric_hist <- function(df, metric, depth, bins = 100) {
  sub <- df %>% filter(num_level_sets == depth)
  ggplot(sub, aes(x = .data[[metric]])) +
    geom_histogram(bins = bins, fill = "#2C7FB8", color = "white", alpha = 0.85) +
    labs(title    = sprintf("%s  |  num_level_sets = %d", metric, depth),
         subtitle = sprintf("n = %d scenarios", nrow(sub)),
         x = metric, y = "count") +
    theme_minimal(base_size = 14)
}


depth_for_filtering <- 3   # choose the CV depth you trust for retention
stats_one_depth <- all_stats_df %>% filter(num_level_sets == depth_for_filtering) %>%
  arrange(scenario_id)



UNCERTAINTY$rmse_logscale <- stats_one_depth$RMSE_log10
UNCERTAINTY$MAPE          <- stats_one_depth$MAPE
UNCERTAINTY$Coverage_50   <- stats_one_depth$Coverage_50

UNCERTAINTY$Coverage_95   <- stats_one_depth$Coverage_95

UNCERTAINTY$r2_logscale <- stats_one_depth$R2_log10
UNCERTAINTY$LogScore_LGOCV <- stats_one_depth$LogScore_LGOCV

UNCERTAINTY$WIS_95   <- stats_one_depth$WIS_95
UNCERTAINTY$WIS_50   <- stats_one_depth$WIS_50

#add the row index before filtering, so each retained scenario carries its 
#original position in UNCERTAINTY (which is what tells you its mesh/formula/priors)
UNCERTAINTY$scenario_id <- seq_len(nrow(UNCERTAINTY))


# examples
plot_metric_hist(all_stats_df, metric = "R2_log10",        depth = 3)
plot_metric_hist(all_stats_df, metric = "RMSE_log10",  depth = 3)
plot_metric_hist(all_stats_df, metric = "MAPE",        depth = 3)
plot_metric_hist(all_stats_df, metric = "Coverage_95", depth = 3)
plot_metric_hist(all_stats_df, metric = "WIS_95", depth = 3)


# define filters by metric
filters <- list(
  RMSE_log10 = UNCERTAINTY %>% filter(rmse_logscale < quantile(rmse_logscale, 0.10)),
  R2_log10       = UNCERTAINTY %>% filter(r2_logscale        > quantile(r2_logscale, 0.90)),
  Coverage95 = UNCERTAINTY %>% filter(Coverage_95 > quantile(Coverage_95, 0.90))
)

#Concretely, build a stacked dataset where each row is a (scenario, retaining-metric) pair:
#keep the original scenario index
labelled <- bind_rows(lapply(names(filters), function(m) {
  filters[[m]] %>% mutate(retain_metric = m)
}))

save(UNCERTAINTY,filters, labelled,mesh_list, formula_list, file="01_Uncertainty_18k_preprocessed.R")
