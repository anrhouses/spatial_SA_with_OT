###############################################################################
# Principal component analysis

#1. Select retained model ensembles based on performance metrics
#2. Build ensemble mean map
#3. PCA explaining more than 70% variance
#4. Accompanying plots

# Author: Corinna Perchtold
###############################################################################

load("01_spatial_points.R")
library(dplyr)
library(ggplot2)
library(sf)
library(scales)
library(patchwork)


out_dir <- "./inla_18000_pred_retained_labelled"

# list all result files
files <- list.files(out_dir, pattern = "^res_scenario_\\d+_.*\\.RData$", full.names = TRUE)
cat("Found", length(files), "result files\n")

# parse scenario id and metric from each filename
meta <- data.frame(
  file        = files,
  scenario_id = as.integer(sub(".*res_scenario_(\\d+)_.*\\.RData$", "\\1", files)),
  retain_metric = sub(".*res_scenario_\\d+_(.*)\\.RData$", "\\1", basename(files)),
  stringsAsFactors = FALSE
)


# build the map matrix: one column per (scenario, metric) result, rows = grid pixels
build_map_matrix_labelled <- function(meta, col_name = "pred_mean") {
  mats <- lapply(meta$file, function(f) {
    e <- new.env()
    load(f, envir = e)
    get(col_name, envir = e)
  })
  M <- do.call(cbind, mats)
  colnames(M) <- paste0(meta$scenario_id, "_", meta$retain_metric)
  M
}

select_metric <- function(meta, metric) {
  sub <- meta[meta$retain_metric == metric, , drop = FALSE]
  if (nrow(sub) == 0L)
    stop("No scenarios for metric '", metric, "'. Available: ",
         paste(unique(meta$retain_metric), collapse = ", "))
  sub
}

meta_rmse <- select_metric(meta, "RMSE_log10")
meta_cov95<-select_metric(meta, "Coverage95")
meta<-rbind(meta_rmse,  meta_cov95)  

# Load and combine
cat("Assembling spatial prediction matrices...")

MAP_PRED_logscale <- build_map_matrix_labelled(meta_rmse, "pred_mean_logscale")
cat("MAP_PRED:", nrow(MAP_PRED_logscale), "pixels x", ncol(MAP_PRED_logscale), "results\n")

cat(" Done.\n")



# ---- coordinates aligned to MAP_PRED rows (pixels) ----
xy <- spatial_points_113AK_df[, c(7, 8)]
colnames(xy) <- c("x", "y")

# ---- collapse the ensemble to ONE map per quantity ----
# ensemble posterior mean map = average of the per-result mean maps
# CV = SD / MEAN * 100
pixel_MEAN <- rowMeans(MAP_PRED_logscale, na.rm = TRUE)

df.plt <- data.frame(x = xy$x, y = xy$y,
                     MEAN = pixel_MEAN)

stations_sf_113AK <- as.data.frame(coords.stations2) %>%
  st_as_sf(coords = c(1, 2), crs = 2154)

# ---- coarsen to a display grid ----
bin_size <- 500   # metres; adjust for resolution of the plotted tiles

df.plt.coarse <- df.plt %>%
  mutate(x_bin = round(x / bin_size) * bin_size,
         y_bin = round(y / bin_size) * bin_size) %>%
  group_by(x_bin, y_bin) %>%
  summarise(MEAN = mean(MEAN, na.rm = TRUE),
            .groups = "drop")

# ---- station + boundary layers ----
st_coords <- st_coordinates(stations_sf_113AK)
df_stations_binned <- data.frame(
  x_bin = round(st_coords[, 1] / bin_size) * bin_size,
  y_bin = round(st_coords[, 2] / bin_size) * bin_size) %>%
  group_by(x_bin, y_bin) %>% summarise(n_stations = n(), .groups = "drop")

station_layer <- geom_point(data = df_stations_binned, aes(x = x_bin, y = y_bin),
                            color = "black", shape = 20, size = 3) #, alpha = 0.8
boundary_layer <- geom_sf(data = entity113AK_boundary, alpha = 0.3,
                          color = "black", linewidth = 0.8)

bold_theme <- theme(
  axis.text = element_text(size = 18, face = "bold"),
  axis.title = element_text(size = 18, face = "bold"),
  legend.text = element_text(size = 18, face = "bold"),
  legend.title = element_text(size = 18, face = "bold"),
  plot.title = element_text(size = 18, face = "bold"),
  plot.subtitle = element_text(size = 18, face = "bold"))

# reusable filler for the warm gradient
warm_fill <- function(v) {
  lo <- min(v, na.rm = TRUE); hi <- max(v, na.rm = TRUE)
  scale_fill_gradientn(
    colors = c("yellow", "orange", "red", "#ca0020"),
    values = scales::rescale(c(lo, hi)),
    limits = c(lo, hi), oob = scales::squish, na.value = "transparent",
    breaks = c( hi*0.25, hi*0.50, hi*0.75, hi),
    labels = round(c( hi*0.25, hi*0.50, hi*0.75, hi), 2))
}

# ---- MEAN map ----
plt.mean <- ggplot() +
  geom_tile(data = df.plt.coarse, aes(x_bin, y_bin, fill = MEAN)) +
  warm_fill(df.plt.coarse$MEAN) +
  station_layer + boundary_layer +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") + 
  scale_x_continuous(
    breaks = c(1.3, 1.8, 2.3, 2.8, 3.3, 3.8),
    labels = c("1.3", "1.8", "2.3", "2.8", "3.3", "3.8")
  )+theme_minimal()+
  scale_y_continuous(
    breaks = c(48.1, 48.6, 49.1, 49.6),
    labels = c("48.1", "48.6", "49.1", "49.6")
  )+
  ggtitle("Posterior mean on log10-scale") + 
  bold_theme

ggsave("pred_mean_logscale.png", plt.mean,
       width = 20, height = 8, dpi = 300, bg = "white")

###################### get PCs explaining more than 70% var ########################

# ============================================================
# FUNCTION: run_pca_pipeline
# Takes a spatial matrix, runs prcomp, and returns:
#   - pca_results : the full prcomp object
#   - BETA        : scenario scores on top PCs
#   - BASIS       : eigenvectors (spatial patterns) of top PCs
#   - n_pcs       : number of PCs needed to reach `var_threshold`
#   - variance_explained / cumulative_variance / labels
#   - M, S        : pixel-wise mean and SD across scenarios
# ============================================================
run_pca_pipeline <- function(MAP_matrix, var_threshold = 0.90) {
  
  # 1. PCA:  scale = TRUE is recommended if columns have different units (e.g. mixing temperature and rainfall)
  pca_results <- prcomp(t(MAP_matrix), center = TRUE, scale = T)
  
  # 2. Variance decomposition
  # This tells you how many "Spatial Patterns" are needed to describe the model differences
  std_dev             <- pca_results$sdev
  variance_explained  <- std_dev^2 / sum(std_dev^2)
  cumulative_variance <- cumsum(variance_explained)
  
  # Find how many PCs explain 90% of the difference between your INLA models
  n_pcs               <- which(cumulative_variance >= var_threshold)[1]
  
  # 3. Extracts the scores (BETA) and the eigenvectors (BASIS) for those top components.
  BETA  <- pca_results$x[, 1:n_pcs, drop = FALSE]          # scenarios × n_pcs
  BASIS <- t(pca_results$rotation[, 1:n_pcs]) # n_pcs × pixels
  
  # 4. Pixel-wise summary statistics
  M <- apply(MAP_matrix, 1, mean)
  S <- apply(MAP_matrix, 1, sd)
  
  list(
    pca_results         = pca_results,
    BETA                = BETA,
    BASIS               = BASIS,
    n_pcs               = n_pcs,
    variance_explained  = variance_explained,
    cumulative_variance = cumulative_variance,
    var_labels          = round(variance_explained  * 100, 2),
    cumvar_labels       = round(cumulative_variance * 100, 2),
    MEAN                = M,
    SD                  = S
  )
}


pca_results_mean_logscale_70 <- run_pca_pipeline(MAP_PRED_logscale, var_threshold = 0.70)

save(pca_results_mean_logscale_70, file="01_PCA_RMSE_log10_Cov95.R")


#For the spatial plot you want to show where in space each PC is dominant — 
#that's the loadings/eigenvectors, which is BASIS:
plot_pca_spatial <- function(pca_results, xy, stations_sf, var_labels,
                             bin_size = 500, white_band = 0.0001, title = "Mean",
                             output_file = "PCA_Spatial_Summary.png") {
  
  n_pcs <- nrow(pca_results$BASIS)
  pc_names <- paste0("PC", 1:n_pcs)
  
  # --- Build fine data frame --------------------------------
  df.plt <- data.frame(x = xy[, 1], y = xy[, 2],
                       MEAN = pca_results$MEAN, 
                       SD   = pca_results$SD,
                       CV = pca_results$SD/pca_results$MEAN*100)
  
  # Add PC columns dynamically
  for (i in 1:n_pcs) {
    df.plt[[paste0("PC", i)]] <- pca_results$BASIS[i, ]
  }
  
  # --- Coarsen to grid --------------------------------------
  df.plt.coarse <- df.plt %>%
    mutate(x_bin = round(x / bin_size) * bin_size,
           y_bin = round(y / bin_size) * bin_size) %>%
    group_by(x_bin, y_bin) %>%
    summarise(MEAN = mean(MEAN, na.rm = TRUE),
              SD   = mean(SD,   na.rm = TRUE),     #SD comes from PCA_result
              CV   = mean(CV, na.rm=T),
              across(all_of(pc_names), mean, na.rm = TRUE),
              .groups = 'drop')
  
  # --- Bin stations -----------------------------------------
  st_coords <- st_coordinates(stations_sf)
  df_stations_binned <- data.frame(
    x_bin = round(st_coords[, 1] / bin_size) * bin_size,
    y_bin = round(st_coords[, 2] / bin_size) * bin_size
  ) %>%
    group_by(x_bin, y_bin) %>%
    summarise(n_stations = n(), .groups = 'drop')
  
  # --- Helper: station layer --------------------------------
  station_layer <- geom_point(
    data  = df_stations_binned, aes(x = x_bin, y = y_bin),
    color = "black", shape = 20, size = 3
  )
  
  boundary_layer <- geom_sf(data = entity113AK_boundary, alpha = 0.3,
                            color = "black", linewidth = 0.8)
  
  
  bold_theme <- theme_minimal() +
    theme(
      axis.text    = element_text(size = 18, face = "bold"),
      axis.title   = element_text(size = 18, face = "bold"),
      legend.text  = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 18, face = "bold"),
      plot.title   = element_text(size = 18, face = "bold"),
      plot.subtitle= element_text(size = 18, face = "bold")
    )
  
  
  mean_min <- min(df.plt.coarse$CV, na.rm = TRUE)
  mean_max <- max(df.plt.coarse$CV, na.rm = TRUE)/6
  
  # --- Mean plot --------------------------------------------
  
  plt.CV <- ggplot() +
    geom_tile(data = df.plt.coarse, aes(x_bin, y_bin, fill = CV)) +
    scale_fill_gradientn(
      colors   = c("yellow", "orange", "red", "#ca0020"),
      values   = scales::rescale(c(0, mean_max/2)),
      limits   = c(mean_min, mean_max),
      oob      = scales::squish,
      na.value = "transparent",
      breaks   = c(mean_min, 
                   mean_max * 0.25, 
                   mean_max * 0.50, 
                   mean_max * 0.75, 
                   mean_max),
      labels   = c(round(mean_min, 0),        # <-- this controls display
                   round(mean_max * 0.25, 0),
                   round(mean_max * 0.50, 0),
                   round(mean_max * 0.75, 0),
                   round(mean_max, 0))
    ) +
    station_layer + boundary_layer +
    scale_x_continuous(
      breaks = c(1.3, 1.8, 2.3, 2.8, 3.3, 3.8),
      labels = c("1.3", "1.8", "2.3", "2.8", "3.3", "3.8")
    )+
    scale_y_continuous(
      breaks = c(48.1, 48.6, 49.1, 49.6),
      labels = c("48.1", "48.6", "49.1", "49.6")
    )+
    labs(x = "Longitude (°E)", y = "Latitude (°N)") +
    ggtitle("Coefficient of variation in %")+bold_theme
  
  mean_min <- min(df.plt.coarse$MEAN, na.rm = TRUE)
  mean_max <- max(df.plt.coarse$MEAN, na.rm = TRUE)
  
  plt.mean <- ggplot() +
    geom_tile(data = df.plt.coarse, aes(x_bin, y_bin, fill = MEAN)) +
    scale_fill_gradientn(
      colors   = c("yellow", "orange", "red", "#ca0020"),
      values   = scales::rescale(c(0, mean_max)),
      limits   = c(mean_min, mean_max),
      oob      = scales::squish,
      na.value = "transparent",
      breaks   = c(mean_min, 
                   mean_max * 0.25, 
                   mean_max * 0.50, 
                   mean_max * 0.75, 
                   mean_max),
      labels   = c(round(mean_min, 0),        # <-- this controls display
                   round(mean_max * 0.25, 0),
                   round(mean_max * 0.50, 0),
                   round(mean_max * 0.75, 0),
                   round(mean_max, 0))
    ) +
    station_layer + boundary_layer +
    scale_x_continuous(
      breaks = c(1.3, 1.8, 2.3, 2.8, 3.3, 3.8),
      labels = c("1.3", "1.8", "2.3", "2.8", "3.3", "3.8")
    )+
    scale_y_continuous(
      breaks = c(48.1, 48.6, 49.1, 49.6),
      labels = c("48.1", "48.6", "49.1", "49.6")
    ) +
    labs(x = "Longitude (°E)", y = "Latitude (°N)")+
    ggtitle(title)+bold_theme
  
  
  
  
  # --- SD plot ----------------------------------------------
  plt.sd <- ggplot() +
    geom_tile(data = df.plt.coarse, aes(x_bin, y_bin, fill = SD)) +
    scale_fill_gradientn(
      colors   = c("yellow", "orange", "red", "#ca0020"),
      values   = scales::rescale(c(0, max(df.plt.coarse$SD, na.rm = TRUE))),
      limits = c(min(df.plt.coarse$SD, na.rm = TRUE), max(df.plt.coarse$SD, na.rm = TRUE)), oob = scales::squish,
      na.value = "transparent") +
    station_layer + boundary_layer +
    scale_x_continuous(
      breaks = c(1.3, 1.8, 2.3, 2.8, 3.3, 3.8),
      labels = c("1.3", "1.8", "2.3", "2.8", "3.3", "3.8")
    )+
    scale_y_continuous(
      breaks = c(48.1, 48.6, 49.1, 49.6),
      labels = c("48.1", "48.6", "49.1", "49.6")
    ) +
    labs(x = "Longitude (°E)", y = "Latitude (°N)")+
    ggtitle("Standard dev.") +bold_theme
  
  # --- PC plots dynamically ---------------------------------
  pc_subtitles <- c("Primary spatial difference across INLA scenarios",
                    "Secondary spatial contrast",
                    "Third spatial contrast",
                    "Fourth spatial contrast")
  
  make_pc_plot <- function(pc_col, pc_num) {
    pc_max <- max(abs(df.plt.coarse[[pc_col]]), na.rm = TRUE)
    white_band <- pc_max * 0.05  # 5% of the max — adjust this percentage
    
    
    ggplot() +
      geom_tile(data = df.plt.coarse, aes(x_bin, y_bin, fill = .data[[pc_col]])) +
      scale_fill_gradient2(
        low      = "#0571b0",   # negative = blue
        mid      = "white",     # zero = white
        high     = "#ca0020",   # positive = red
        midpoint = 0.00,           # white is exactly at 0
        #limits   = c(-pc_max, pc_max),
        oob      = scales::squish,
        na.value = "transparent"
      )+
      station_layer + boundary_layer + scale_x_continuous(
        breaks = c(1.3, 1.8, 2.3, 2.8, 3.3, 3.8),
        labels = c("1.3", "1.8", "2.3", "2.8", "3.3", "3.8")
      )+
      scale_y_continuous(
        breaks = c(48.1, 48.6, 49.1, 49.6),
        labels = c("48.1", "48.6", "49.1", "49.6")
      ) +
      labs(x = "Longitude (°E)", y = "Latitude (°N)")+
      labs(title    = paste0("PC", pc_num, " (", var_labels[pc_num], "% variance)"),
           subtitle = pc_subtitles[pc_num])  +bold_theme
  }
  
  pc_plots <- lapply(1:n_pcs, function(i) make_pc_plot(paste0("PC", i), i))
  
  # --- Save: mean + sd + all PCs ----------------------------
  all_plots  <- c(plt.CV, pc_plots)#list(plt.mean, plt.sd), pc_plots)
  n_plots    <- length(all_plots)
  n_cols     <- 2
  n_rows     <- ceiling(n_plots / n_cols)
  
  png(output_file, width = 2600 * n_cols, height = 2200 * n_rows, res = 400)
  gridExtra::grid.arrange(grobs = all_plots, ncol = n_cols)
  dev.off()
  
  message("Saved: ", output_file, " (", n_pcs, " PCs, ", n_plots, " panels)")
}

stations_sf_113AK <- as.data.frame(coords.stations2) %>%
  st_as_sf(coords = c(1, 2), crs = 2154)

plot_pca_spatial(pca_results_mean_logscale_70,  xy, stations_sf_113AK, pca_results_mean_logscale_70$var_labels,title="Mean Prediction on Log-Scale", output_file = "PCA_mean_logscale_113AK_with_scale.png")

