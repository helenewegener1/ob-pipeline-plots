#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(viridis)
  library(scales) 
  library(jsonlite)
  library(stringr)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--perf_input"), type = "character", default = NULL, 
              help = "Path to the performances.tsv file", metavar = "character"),
  make_option(c("-j", "--meta_json"), type = "character", default = NULL, 
              help = "Path to dataset_metadata.json", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "Figure4_Runtime_Metadata.png", 
              help = "Output filename", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validation
if (is.null(opt$perf_input) || is.null(opt$meta_json)) {
  print_help(opt_parser)
  stop("Both --perf_input and --meta_json must be provided.", call. = FALSE)
}

path_perf <- opt$perf_input
path_meta_json <- opt$meta_json
output_file <- opt$output

# ------------------------------------------------------------------------------
# 1. DATA LOADING
# ------------------------------------------------------------------------------
df_perf_raw <- read_tsv(path_perf, show_col_types = FALSE)
metadata_json <- fromJSON(path_meta_json)

# ------------------------------------------------------------------------------
# 2. MAPPINGS & REFERENCE DATA
# ------------------------------------------------------------------------------
name_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 'CV',
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 'HT',
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 'COVID',
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 'SB',
  'dataset_name-Samusik_seed-42'                    = 'MBM',
  'dataset_name-Transformed_seed-42'                = 'TF',
  'dataset_name-flowcyt_seed-42'                    = 'HBM',
  'dataset_name-Levine_seed-42'                     = 'LV',
  "dataset_name-FR-FCM-ZZRQ_seed-42"                = "DCI"
)

model_map <- c(
  'cyanno'      = "CyAnno",
  'cygate'      = "CyGATE",
  'dgcytof'     = "DGCytof",
  'gatemeclass' = "GateMeClass",
  'lda'         = "CyTOF LC",
  'random'      = "Random"
)

# Robust Marker Map
marker_map_df <- data.frame(
  dataset_id = c(
    'dataset_name-FR-FCM-Z238_infection_final_seed-42',
    'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42',
    'dataset_name-FR-FCM-Z2KP_virus_final_seed-42',
    'dataset_name-FR-FCM-Z3YR_seed-42',
    'dataset_name-Samusik_seed-42',
    'dataset_name-Transformed_seed-42',
    'dataset_name-flowcyt_seed-42',
    'dataset_name-Levine_seed-42',
    "dataset_name-FR-FCM-ZZRQ_seed-42"
  ),
  n_markers = c(37, 24, 24, 38, 39, 33, 12, 32, 9),
  stringsAsFactors = FALSE
)

# --- COLORS (Set1 Palette - Consistent with Fig 3) ---
tool_colors <- c(
  "CyAnno"      = "#E41A1C",  # Red
  "CyGATE"      = "#377EB8",  # Blue
  "DGCytof"     = "#4DAF4A",  # Green
  "GateMeClass" = "#984EA3",  # Purple
  "CyTOF LC"    = "#FF7F00",  # Orange
  "Random"      = "#525252"   # Dark Grey
)

# ------------------------------------------------------------------------------
# 3. DATA PROCESSING
# ------------------------------------------------------------------------------

# A. Parse Performance Data & Average Runtime
# Note: Adapting regex to handle dataset naming conventions in params
df_time_avg <- df_perf_raw %>%
  # Exclude internal pipeline steps
  filter(!module %in% c("data_import", "data_preprocessing", "flow_metrics", "f1_score", "metric_collector", "metrics")) %>%
  mutate(
    # Extract dataset name
    extracted_name = str_match(params, 'dataset_name[^:]+:\\s*\"+([^\"]+)')[,2],
    # Extract seed
    extracted_seed = str_match(params, 'seed[^:]+:\\s*\"+([^\"]+)')[,2]
  ) %>%
  # Construct ID to match metadata keys
  mutate(dataset_id = paste0("dataset_name-", extracted_name, "_seed-", extracted_seed)) %>%
  filter(!is.na(extracted_name)) %>%
  # Apply Model Map
  mutate(model = recode(module, !!!model_map)) %>%
  # Average runtime across seeds/CVs
  group_by(dataset_id, model) %>%
  summarise(mean_time_sec = mean(s, na.rm = TRUE), .groups = "drop")

# B. Metadata Extraction Function
extract_metadata_stats <- function(meta_list) {
  ds_ids <- names(meta_list)
  results <- lapply(ds_ids, function(id) {
    entry <- meta_list[[id]]
    cells <- as.numeric(unlist(entry$cells_per_sample))
    avg_cells <- if(all(is.na(cells))) NA else mean(cells, na.rm = TRUE)
    
    data.frame(
      dataset_id    = id,
      n_samples     = as.numeric(entry$sample_count),
      n_populations = as.numeric(entry$population_count),
      mean_cells    = avg_cells,
      stringsAsFactors = FALSE
    )
  })
  return(do.call(rbind, results))
}

df_metadata <- extract_metadata_stats(metadata_json)

# C. Merge and Calculate
df_plot <- df_time_avg %>%
  # 1. Join Metadata
  left_join(df_metadata, by = "dataset_id") %>%
  # 2. Join Markers
  left_join(marker_map_df, by = "dataset_id") %>%
  # 3. Filter
  filter(!str_detect(dataset_id, regex("sub-sampling", ignore_case = TRUE))) %>%
  filter(!str_detect(dataset_id, regex("Levine", ignore_case = TRUE))) %>%
  filter(!str_detect(model, regex("random", ignore_case = TRUE))) %>%
  # 4. Clean Names
  mutate(dataset_clean = recode(dataset_id, !!!name_map)) %>%
  filter(!is.na(dataset_clean)) %>%
  # --- NEW CALCULATION: Markers per Population ---
  mutate(markers_per_pop = n_markers / n_populations)

# ------------------------------------------------------------------------------
# 4. PLOTTING
# ------------------------------------------------------------------------------

# Shared Theme
theme_gb_scatter <- theme_bw(base_size = 9) +
  theme(
    text = element_text(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    legend.position = "none"
  )

create_time_plot <- function(data, x_var, x_label, is_log_x = FALSE) {
  
  # Ensure data is sorted by x_var so lines connect in order
  data <- data %>% arrange(.data[[x_var]])
  
  p <- ggplot(data, aes(x = .data[[x_var]], y = mean_time_sec, color = model, group = model)) +
    
    # Connect dots (geom_line) instead of regression
    geom_line(linewidth = 0.5, alpha = 0.5) +
    
    # Scatter Points
    geom_point(aes(fill = model), shape = 21, color = "black", stroke = 0.2, size = 2.5, alpha = 0.8) +
    
    # --- COLORS ---
    scale_color_manual(values = tool_colors, name = "Method") +
    scale_fill_manual(values = tool_colors, name = "Method") +
    
    # --- Y AXIS: LOG SCALE FOR RUNTIME ---
    scale_y_log10(
      labels = label_log(), 
      breaks = trans_breaks("log10", function(x) 10^x)
    ) +
    
    labs(x = x_label, y = "Mean Runtime (s)") +
    theme_gb_scatter
  
  if(is_log_x) {
    p <- p + scale_x_log10(
      labels = label_log(), 
      breaks = trans_breaks("log10", function(x) 10^x)
    )
  }
  return(p)
}

# Generate Specific Plots

# Plot 1: Mean Cells (Log X) vs Runtime (Log Y)
p_cells <- create_time_plot(df_plot, "mean_cells", "Mean Cells / Sample", is_log_x = TRUE)

# Plot 2: Markers per Pop (Linear X) vs Runtime (Log Y)
p_markers_pop <- create_time_plot(df_plot, "markers_per_pop", "Number of Markers / Population")

# ------------------------------------------------------------------------------
# 5. ASSEMBLE AND SAVE
# ------------------------------------------------------------------------------

# Side by side layout
final_fig <- (p_cells + p_markers_pop) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a') & 
  theme(
    legend.position = "bottom", 
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 9),
    legend.key.size = unit(4, "mm"),
    plot.tag = element_text(face = "bold", size = 12)
  )

# Ensure output directory exists (basic check)
if(!str_detect(output_file, "\\.png$")) {
  output_file <- paste0(output_file, ".png")
}

# Reduced height for single row
ggsave(output_file, 
       plot = final_fig, 
       width = 180, 
       height = 90, 
       units = "mm", 
       dpi = 600)

print(paste("Saved Figure 4 to", output_file))

###
# Usage example:
#Rscript Fig4_plot.r \
#  --perf_input ../ob-blob-metrics/out/performances.tsv \
#  --meta_json ../ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
#  --output ../ob-pipeline-plots/Figure4_Runtime_Metadata.png
###