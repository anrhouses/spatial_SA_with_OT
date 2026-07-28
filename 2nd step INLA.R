###############################################################################
# Model fitting and prediction function

# For each retained scenario:
#1. Rebuild the SPDE model
#2. Fit against all observations
#3. Predict over the complete prediction domain
#4. Save posterior summaries

# Author: Corinna Perchtold
###############################################################################

library(ggplot2)
library(sf)
library(sp)
library(terra)
library(fmesher)
library(Matrix)
library(gstat)

library(INLA)
library(dplyr)

load("01_Uncertainty_18k_preprocessed.R")
set.seed(1234)

# integer code for the metric
labelled$metric_code <- match(labelled$retain_metric, names(filters))

# the ORIGINAL scenario row-index in UNCERTAINTY 
stopifnot("scenario_id" %in% names(labelled))

cat(sprintf("labelled rows (scenario x metric): %d  (unique scenarios: %d)\n",
            nrow(labelled), length(unique(labelled$scenario_id))))

out_dir <- "./inla_18000_pred_retained_labelled"
if (!dir.exists(out_dir)) dir.create(out_dir)



run_inla_labelled <- function(j, LAB) {
  row    <- LAB[j, ]
  sid    <- row$scenario_id
  metric <- row$retain_metric
  
  # checkpoint tagged by scenario AND metric so duplicates don't collide
  file_path <- file.path(out_dir, sprintf("res_scenario_%05d_%s.RData", sid, metric))
  if (file.exists(file_path)) return(NULL)
  
  curr_mesh <- mesh_list[[row$mesh_idx]]
  curr_form <- formula_list[[row$form_idx]]
  
  spde <- inla.spde2.pcmatern(
    mesh        = curr_mesh,
    prior.range = c(row$prior_range, 0.5),
    prior.sigma = c(row$prior_sigma, 0.5),  
    alpha       = row$alpha_value)
  
  A_obs  <- inla.spde.make.A(curr_mesh, loc = coords.stations2)
  A_pred <- inla.spde.make.A(curr_mesh, loc = coords.prediction2)
  index_space <- inla.spde.make.index("spatial_field", n.spde = spde$n.spde)
  
  stk_obs <- inla.stack(
    data = list(y = log10(obs.data2$resultat)),
    A = list(A_obs, 1),
    effects = list(c(index_space, list(Intercept = 1)), list(covars.data2)),
    tag = "obs")
  
  stk_pred <- inla.stack(
    data = list(y = NA),
    A = list(A_pred, 1),
    effects = list(c(index_space, list(Intercept = 1)),
                   list(spatial_points_113AK_df[, c(1,2,3,4,5,6)])),
    tag = "pred")
  
  join <- inla.stack(stk_obs, stk_pred)
  
  res <- inla(curr_form, family = "gaussian",
              data = inla.stack.data(join, spde = spde),
              quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
              control.predictor = list(A = inla.stack.A(join), compute = TRUE, link = 1),
              control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = FALSE,
                                     return.marginals.predictor = FALSE),
              inla.mode = "classic",
              control.inla = list(int.strategy = "eb", strategy = "adaptive"),
              safe = TRUE, verbose = F)
  
  pidx <- inla.stack.index(join, "pred")$data
  fv   <- res$summary.fitted.values
  
  pred_mean_logscale     <- fv[pidx, "mean"]
  pred_sd_logscale       <- fv[pidx, "sd"]
  pred_quant025_logscale <- fv[pidx, "0.025quant"]
  pred_quant25_logscale  <- fv[pidx, "0.25quant"]
  pred_quant50_logscale  <- fv[pidx, "0.5quant"]
  pred_quant75_logscale  <- fv[pidx, "0.75quant"]
  pred_quant975_logscale <- fv[pidx, "0.975quant"]
  
  # back-transforms (lognormal mean + lognormal sd)
  pred_mean <- 10^(pred_mean_logscale + (pred_sd_logscale^2 * log(10)) / 2)
  pred_sd   <- sqrt((10^(pred_sd_logscale^2 * log(10)) - 1) *
                      10^(2 * pred_mean_logscale + pred_sd_logscale^2 * log(10)))
  pred_quant025 <- 10^(pred_quant025_logscale)
  pred_quant25  <- 10^(pred_quant25_logscale)
  pred_quant50  <- 10^(pred_quant50_logscale)
  pred_quant75  <- 10^(pred_quant75_logscale)
  pred_quant975 <- 10^(pred_quant975_logscale)
  
  save(pred_mean, pred_sd, pred_quant025, pred_quant25, pred_quant50,
       pred_quant75, pred_quant975,
       pred_mean_logscale, pred_sd_logscale, pred_quant025_logscale,
       pred_quant25_logscale, pred_quant50_logscale,
       pred_quant75_logscale, pred_quant975_logscale,
       file = file_path)
  
  rm(res, join); gc()
  return(NULL)   # results are on disk; no need to ship back
}

library(tictoc); library(future); library(future.apply); library(progressr)
handlers(handler_txtprogressbar(style = 3))

total_cores <- parallel::detectCores(logical = FALSE)
inla_threads <-as.integer(Sys.getenv("INLA_NUM_THREADS", unset = "2"))
n_workers <- floor(total_cores / inla_threads) 
cat(sprintf("Parallelisation: %d physical cores, %d INLA threads -> %d parallel workers\n",
            total_cores, inla_threads, n_workers))


# Force each INLA fit to respect the thread budget (prevents 64*128 explosion)
INLA::inla.setOption(num.threads = paste0(inla_threads, ":1"))

# Cap future's per-worker object transfer so it errors loudly instead of silently
options(future.globals.maxSize = 8 * 1024^3)   # 8 GB ceiling per worker

plan(multisession, workers = n_workers)

n_rows <- nrow(labelled)
cat(sprintf("Predicting %d (scenario x metric) rows\n", n_rows))

tic()
with_progress({
  future_lapply(1:n_rows, function(j) {
    library(INLA)
    tryCatch(run_inla_labelled(j, labelled),
             error = function(e)
               write(sprintf("row %d (scenario %d, %s) failed: %s",
                             j, labelled$scenario_id[j], labelled$retain_metric[j], e$message),
                     file = "pred_error_labelled.txt", append = TRUE))
    
  }, future.seed = TRUE)
})
toc()
plan(sequential)
cat("Prediction done.\n")