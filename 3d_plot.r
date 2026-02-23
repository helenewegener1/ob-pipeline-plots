#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 1. LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(stringr)
  library(scatterplot3d)
})

# ------------------------------------------------------------------------------
# 2. COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-f", "--f1_macro"), type = "character", default = NULL,
              help = "Path to f1_macro_by_crossvalidation.tsv", metavar = "character"),
  make_option(c("-p", "--perf_input"), type = "character", default = NULL,
              help = "Path to performances.tsv", metavar = "character"),
  make_option(c("-j", "--meta_json"), type = "character", default = NULL,
              help = "Path to dataset_metadata.json", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "Figure_3D_Performance.png",
              help = "Output filename (.png)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$f1_macro) || is.null(opt$perf_input) || is.null(opt$meta_json)) {
  print_help(opt_parser)
  stop("All inputs (--f1_macro, --perf_input, --meta_json) must be provided.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 3. GLOBAL MAPPINGS & COLORS
# ------------------------------------------------------------------------------
name_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 'CV',
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 'HT',
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 'COVID',
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 'SB',
  'dataset_name-Samusik_seed-42'                    = 'MBM',
  'dataset_name-Transformed_seed-42'                = 'TF',
  'dataset_name-flowcyt_seed-42'                    = 'HBM',
  'dataset_name-Levine_seed-42'                     = 'LV'
)

populations_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 16,
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 6,
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 6,
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 28,
  'dataset_name-Samusik_seed-42'                    = 15,
  'dataset_name-Transformed_seed-42'                = 8,
  'dataset_name-flowcyt_seed-42'                    = 5,
  'dataset_name-Levine_seed-42'                     = 15
)

marker_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 37,
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 24,
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 24,
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 38,
  'dataset_name-Samusik_seed-42'                    = 39,
  'dataset_name-Transformed_seed-42'                = 33,
  'dataset_name-flowcyt_seed-42'                    = 12,
  'dataset_name-Levine_seed-42'                     = 32
)

model_map <- c(
  'cyanno' = "CyAnno",
  'cygate' = "CyGATE",
  'dgcytof' = "DGCyTOF",
  'gatemeclass[V]' = "GateMeClass-V",
  'gatemeclass[E]' = "GateMeClass-E",
  'lda' = "CyTOF LC",
  'knn' = "KNN",
  'random' = "Random"
)

tool_colors <- c(
  "CyAnno"        = "#E41A1C",  
  "CyGATE"        = "#377EB8",  
  "DGCyTOF"       = "#4DAF4A",  
  "CyTOF LC"      = "#FF7F00",  
  "GateMeClass-E" = "#6A3D9A",  
  "GateMeClass-V" = "#CAB2D6",  
  "KNN"           = "#17BECF",  
  "Random"        = "#525252"   
)

# ------------------------------------------------------------------------------
# 4. DATA PROCESSING
# ------------------------------------------------------------------------------
message("Loading data...")
df_f1_raw     <- read_tsv(opt$f1_macro, show_col_types = FALSE)
df_perf_raw   <- read_tsv(opt$perf_input, show_col_types = FALSE)
metadata_json <- fromJSON(opt$meta_json)

# Setup Metadata
cell_stats <- do.call(rbind, lapply(names(metadata_json), function(id) {
  cells <- as.numeric(unlist(metadata_json[[id]]$cells_per_sample))
  data.frame(dataset_id = id, mean_cells = ifelse(all(is.na(cells)), NA, mean(cells, na.rm = TRUE)))
}))

df_general_meta <- data.frame(dataset_id = names(name_map)) %>%
  mutate(n_populations = recode(dataset_id, !!!populations_map),
         n_markers     = recode(dataset_id, !!!marker_map)) %>%
  left_join(cell_stats, by = "dataset_id") %>%
  mutate(markers_per_pop = n_markers / n_populations)

# Process F1 Macro
df_f1_avg <- df_f1_raw %>%
  filter(!is.na(f1_macro), f1_macro > 0) %>%
  filter(!str_detect(dataset, "panel_CD20|panel_CD56")) %>%
  mutate(model = recode(model, !!!model_map)) %>%
  filter(model != "Random") %>%
  group_by(dataset, model) %>%
  summarise(mean_f1_macro = mean(f1_macro, na.rm = TRUE), .groups = "drop")

# Process Performance
df_time_avg <- df_perf_raw %>%
  filter(!module %in% c("data_import", "data_preprocessing", "flow_metrics", "f1_score", "metric_collector", "metrics")) %>%
  mutate(
    extracted_name = str_match(params, 'dataset_name[^:]+:\\s*\"+([^\"]+)')[, 2],
    extracted_seed = str_match(params, 'seed[^:]+:\\s*\"+([^\"]+)')[, 2],
    dataset_id = paste0("dataset_name-", extracted_name, "_seed-", extracted_seed),
    model = recode(module, !!!model_map)
  ) %>%
  filter(!is.na(extracted_name), model != "Random") %>%
  # Sync Perf with F1 (Drops runs that failed F1)
  semi_join(df_f1_avg, by = c("dataset_id" = "dataset", "model" = "model")) %>%
  group_by(dataset_id, model) %>%
  summarise(mean_time_sec = mean(s, na.rm = TRUE), .groups = "drop")

# Final Plot Tables
df_plot_f1 <- df_f1_avg %>% 
  inner_join(df_general_meta, by = c("dataset" = "dataset_id")) %>%
  mutate(log_cells = log10(mean_cells), color = tool_colors[model])

df_plot_time <- df_time_avg %>% 
  inner_join(df_general_meta, by = "dataset_id") %>%
  mutate(log_cells = log10(mean_cells), log_time = log10(mean_time_sec), color = tool_colors[model])

# Filter tool_colors to only those present in the final dataset for the legend
active_models <- sort(unique(df_plot_f1$model))
active_colors <- tool_colors[active_models]

# ------------------------------------------------------------------------------
# 5. 3D PLOTTING & SAVING
# ------------------------------------------------------------------------------
message("Generating 3D plot...")

png(filename = opt$output, width = 14, height = 6.5, units = "in", res = 300)

layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
par(oma = c(2, 2, 3, 10)) # Margins to leave room for the legend on the right

# --- PANEL A: F1 Macro Score ---
par(mar = c(4, 3, 2, 2))
s3d_f1 <- scatterplot3d(
  x = df_plot_f1$log_cells, y = df_plot_f1$markers_per_pop, z = df_plot_f1$mean_f1_macro,
  color = df_plot_f1$color,
  pch = 16, type = "h", lty.hplot = 3, angle = 45, scale.y = 0.8,
  zlim = c(0, 1), 
  main = "Model Accuracy Profile",
  xlab = "Log10(Mean Cells)", ylab = "Markers / Population", zlab = "Mean Macro F1"
)

# --- PANEL B: Runtime ---
par(mar = c(4, 3, 2, 2))
s3d_time <- scatterplot3d(
  x = df_plot_time$log_cells, y = df_plot_time$markers_per_pop, z = df_plot_time$log_time,
  color = df_plot_time$color,
  pch = 16, type = "h", lty.hplot = 3, angle = 45, scale.y = 0.8,
  main = "Model Efficiency Profile",
  xlab = "Log10(Mean Cells)", ylab = "Markers / Population", zlab = "Log10(Runtime in seconds)"
)

# --- UNIFIED LEGEND ---
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(
  "right", 
  legend = active_models, 
  col = active_colors, 
  pch = 16, 
  bty = "n", 
  cex = 1.2,
  title = "Method",
  inset = c(0.02, 0)
)

dev.off()
message(paste("Saved 3D combined figure to:", opt$output))

### Example usage:
#Rscript ./3d_plot.r \
#  --f1_macro ../ob-blob-metrics/out/metric_collectors/metrics_report/f1_macro_by_crossvalidation.tsv \
#  --perf_input ../ob-blob-metrics/out/performances.tsv \
#  --meta_json ../ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
#  --output ../ob-pipeline-plots/Figure3_3D_Performance.png