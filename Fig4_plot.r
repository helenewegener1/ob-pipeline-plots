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
# 2. METADATA EXTRACTION
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
# ------------------------------------------------------------------------------
# 3. PARSE PERFORMANCE & SUMMARIZE
# ------------------------------------------------------------------------------
# A. Parse Performance Data & Average Runtime
df_time_avg <- df_perf_raw %>%
  # Exclude internal pipeline steps
  filter(!module %in% c("data_import", "data_preprocessing", "flow_metrics", "f1_score", "metric_collector", "metrics")) %>%
  # Regex to extract dataset name and seed from the params string
  mutate(
    extracted_name = str_match(params, 'dataset_name[^:]+:\\s*\"+([^\"]+)')[,2],
    extracted_seed = str_match(params, 'seed[^:]+:\\s*\"+([^\"]+)')[,2],
    dataset_id = paste0("dataset_name-", extracted_name, "_seed-", extracted_seed)
  ) %>%
  filter(!is.na(extracted_name)) %>%
  # Apply Model Map immediately
  mutate(model = recode(module, !!!model_map)) %>%
  # Average runtime across seeds/CVs
  group_by(dataset_id, model) %>%
  summarise(mean_time_sec = mean(s, na.rm = TRUE), .groups = "drop")

df_metadata <- extract_metadata_stats(metadata_json)

# C. Merge
df_plot <- df_time_avg %>%
  left_join(df_metadata, by = "dataset_id") %>%
  left_join(marker_map_df, by = "dataset_id") %>%
  filter(!str_detect(dataset_id, regex("sub-sampling", ignore_case = TRUE))) %>%
  filter(!str_detect(dataset_id, regex("Levine", ignore_case = TRUE))) %>%
  filter(!str_detect(model, regex("random", ignore_case = TRUE))) %>%
  mutate(dataset_clean = recode(dataset_id, !!!name_map)) %>%
  filter(!is.na(dataset_clean))

# 3. THEME & PLOTTING FUNCTION

# Shared Theme (Identical to Fig 3)
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
  
  p <- ggplot(data, aes(x = .data[[x_var]], y = mean_time_sec, color = model, fill = model)) +
    # Linear Trend Line
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, alpha = 0.5) +
    
    # Scatter Points
    geom_point(shape = 21, color = "black", stroke = 0.2, size = 2.5, alpha = 0.8) +
    
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

# 4. GENERATE PLOTS

p_cells   <- create_time_plot(df_plot, "mean_cells", "Mean Cells / Sample", is_log_x = TRUE)
p_markers <- create_time_plot(df_plot, "n_markers", "Number of Markers")
p_samples <- create_time_plot(df_plot, "n_samples", "Number of Samples")
p_pops    <- create_time_plot(df_plot, "n_populations", "Number of Populations")

# 5. ASSEMBLE (Unified Legend)
final_fig <- (p_cells + p_markers) / (p_samples + p_pops) + 
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

# 6. SAVE
ggsave(output_file, 
       plot = final_fig, 
       width = 180, 
       height = 160, 
       units = "mm", 
       dpi = 600)

print(paste("Saved Figure 4 to", output_file))

###
# Usage example:
#Rscript Fig4_plot.r \
#  --perf_input ./ob-blob-metrics/out/performances.tsv \
#  --meta_json ./ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
#  --output ./ob-pipeline-plots/Figure4_Runtime_Metadata.png
###