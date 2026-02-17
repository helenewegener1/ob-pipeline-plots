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
  make_option(c("-i", "--f1_weighted"), type = "character", default = NULL, 
              help = "Path to samples_vs_f1_weighted.tsv", metavar = "character"),
  make_option(c("-j", "--meta_json"), type = "character", default = NULL, 
              help = "Path to dataset_metadata.json", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "Figure3_Metadata_Performance.png", 
              help = "Output filename", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validation
if (is.null(opt$f1_weighted) || is.null(opt$meta_json)) {
  print_help(opt_parser)
  stop("Both --f1_weighted and --meta_json must be provided.", call. = FALSE)
}

path_f1 <- opt$f1_weighted
path_meta_json <- opt$meta_json
output_file <- opt$output

# ------------------------------------------------------------------------------
# 1. DATA LOADING & MAPPINGS
# ------------------------------------------------------------------------------
df_f1_raw <- read_tsv(path_f1, show_col_types = FALSE)
metadata_json <- fromJSON(path_meta_json)

name_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 'CV',
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 'HT',
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 'COVID',
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 'SB',
  'dataset_name-Samusik_seed-42'                    = 'MBM',
  'dataset_name-Transformed_seed-42'                = 'TF',
  'dataset_name-flowcyt_seed-42' = 'HBM',
  'dataset_name-Levine_seed-42' = 'LV',
  "dataset_name-panel_CD20_seed-42" = "DCI-CD20",
  "dataset_name-panel_CD56_seed-42" = "DCI-CD56"
)

model_map <- c(
'cyanno' = "CyAnno",
'cygate' = "CyGATE",
'dgcytof' = "DGCytof",
'gatemeclass' = "GateMeClass",
'lda' = "CyTOF LC",
'knn' = "KNN",
'random' = "Random"
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
    "dataset_name-panel_CD20_seed-42",
    "dataset_name-panel_CD56_seed-42"
  ),
  n_markers = c(37, 24, 24, 38, 39, 33, 12, 32, 8,8),
  stringsAsFactors = FALSE
)

# --- COLORS (Set1 Palette) ---
tool_colors <- c(
  "CyAnno"      = "#E41A1C",  # Bold Red
  "CyGATE"      = "#377EB8",  # Strong Blue
  "DGCytof"     = "#4DAF4A",  # Vivid Green
  "GateMeClass" = "#984EA3",  # Deep Purple
  "CyTOF LC"    = "#FF7F00",  # Strong Orange
  'KNN' = "#ec7ed0",
  "Random"      = "#525252"   # Dark Grey (Charcoal)
)

# ------------------------------------------------------------------------------
# 2. PROCESSING
# ------------------------------------------------------------------------------
# A. Summarize F1 Weighted Median per (Dataset x Model)
df_f1_avg <- df_f1_raw %>%
  mutate(model = recode(model, !!!model_map)) %>% 
  group_by(dataset, model) %>%
  summarise(mean_f1_weighted = mean(f1_weighted_median, na.rm = TRUE), .groups = "drop")

# B. Metadata Extraction Function
extract_metadata_stats <- function(meta_list) {
  ds_ids <- names(meta_list)
  
  results <- lapply(ds_ids, function(id) {
    entry <- meta_list[[id]]
    
    # Calculate mean cells per sample
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

# C. Merge and Final Cleanup
df_plot <- df_f1_avg %>%
  left_join(df_metadata, by = c("dataset" = "dataset_id")) %>%
  left_join(marker_map_df, by = c("dataset" = "dataset_id")) %>%
  filter(
    !str_detect(dataset, regex("sub-sampling|Levine", ignore_case = TRUE)),
    !str_detect(model, regex("random", ignore_case = TRUE))
  ) %>%
  mutate(
    dataset_clean = recode(dataset, !!!name_map),
    markers_per_pop = n_markers / n_populations
  ) %>%
  filter(!is.na(dataset_clean),mean_f1_weighted > 0) %>%
  # drop_na is a concise way to handle multiple columns at once
  tidyr::drop_na(mean_f1_weighted, mean_cells, markers_per_pop)

# 3. THEME & PLOTTING FUNCTION

# Genome Biology Clean Theme
theme_gb_scatter <- theme_bw(base_size = 9) +
  theme(
    text = element_text(color = "black"),
    # Clean Borders
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    
    # Axis Text
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    
    # Legend (Hidden in individual plots, collected later)
    legend.position = "none"
  )

create_plot <- function(data, x_var, x_label, is_log = FALSE) {
  
  p <- ggplot(data, aes(x = .data[[x_var]], y = mean_f1_weighted, color = model, group = model)) +
    
    # Connect dots (Instead of regression)
    geom_line(linewidth = 0.5, alpha = 0.5) +
    
    # Scatter Points
    geom_point(aes(fill = model), shape = 21, color = "black", stroke = 0.2, size = 2.5, alpha = 0.8) +
    
    # --- APPLY CONSISTENT COLORS ---
    scale_color_manual(values = tool_colors, name = "Method") +
    scale_fill_manual(values = tool_colors, name = "Method") +
    
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = x_label, y = "Mean F1-Score") +
    theme_gb_scatter
  
  if(is_log) {
    # Log scale for cells/sample
    p <- p + scale_x_log10(
      labels = label_log(), 
      breaks = trans_breaks("log10", function(x) 10^x)
    )
  }
  return(p)
}

# 4. GENERATE PLOTS

# Plot 1: Mean Cells per Sample (Log Scale)
p_cells <- create_plot(df_plot, "mean_cells", "Mean Cells / Sample", is_log = TRUE)

# Plot 2: Markers per Population (New Request)
p_markers_pop <- create_plot(df_plot, "markers_per_pop", "Number of Markers / Population")

# 5. ASSEMBLE WITH PATCHWORK
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

# 6. SAVE
# Ensure output directory exists if filename implies one, otherwise standard save
outfile <- output_file
if(!str_detect(outfile, "\\.png$")) {
  outfile <- paste0(outfile, ".png")
}

# Adjusted height since we only have 1 row of plots now
ggsave(outfile, 
       plot = final_fig, 
       width = 180, 
       height = 90, # Reduced height for single row
       units = "mm", 
       dpi = 600)

print(paste("Saved Figure 3 to", outfile))

###
# Usage example:
#Rscript Fig3_plot.r \
#  --f1_weighted ../ob-blob-metrics/out/metric_collectors/metrics_report/samples_vs_f1_weighted.tsv \
#  --meta_json ../ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
#  --output ../ob-pipeline-plots/Figure3_Performance.png
###