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
  library(cowplot)   # <- instead of patchwork
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-f", "--f1_weighted"), type = "character", default = NULL,
              help = "Path to samples_vs_f1_weighted.tsv", metavar = "character"),
  make_option(c("-p", "--perf_input"), type = "character", default = NULL,
              help = "Path to performances.tsv", metavar = "character"),
  make_option(c("-j", "--meta_json"), type = "character", default = NULL,
              help = "Path to dataset_metadata.json", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "Figure_Combined_Grid.png",
              help = "Output filename", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$f1_weighted) || is.null(opt$perf_input) || is.null(opt$meta_json)) {
  print_help(opt_parser)
  stop("All inputs (--f1_weighted, --perf_input, --meta_json) must be provided.", call. = FALSE)
}

path_f1 <- opt$f1_weighted
path_perf <- opt$perf_input
path_meta_json <- opt$meta_json
output_file <- opt$output

# ------------------------------------------------------------------------------
# 1. GLOBAL MAPPINGS & THEMES
# ------------------------------------------------------------------------------

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
  'KNN' = "#ec7ed0"
)

theme_gb_scatter <- theme_bw(base_size = 9) +
  theme(
    text = element_text(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black")
  )

# ------------------------------------------------------------------------------
# 2. DATA LOADING & PRE-PROCESSING
# ------------------------------------------------------------------------------
df_f1_raw <- read_tsv(path_f1, show_col_types = FALSE)
df_perf_raw <- read_tsv(path_perf, show_col_types = FALSE)
metadata_json <- fromJSON(path_meta_json)

extract_metadata_stats <- function(meta_list) {
  ds_ids <- names(meta_list)
  results <- lapply(ds_ids, function(id) {
    entry <- meta_list[[id]]
    cells <- as.numeric(unlist(entry$cells_per_sample))
    avg_cells <- if (all(is.na(cells))) NA else mean(cells, na.rm = TRUE)

    data.frame(
      dataset_id    = id,
      n_samples     = as.numeric(entry$sample_count),
      n_populations = as.numeric(entry$population_count),
      mean_cells    = avg_cells,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}
df_metadata <- extract_metadata_stats(metadata_json)

df_f1_avg <- df_f1_raw %>%
  mutate(model = recode(model, !!!model_map)) %>%
  group_by(dataset, model) %>%
  summarise(mean_f1_weighted = mean(f1_weighted_median, na.rm = TRUE), .groups = "drop")

df_time_avg <- df_perf_raw %>%
  filter(!module %in% c("data_import", "data_preprocessing", "flow_metrics", "f1_score",
                       "metric_collector", "metrics")) %>%
  mutate(
    extracted_name = str_match(params, 'dataset_name[^:]+:\\s*\"+([^\"]+)')[, 2],
    extracted_seed = str_match(params, 'seed[^:]+:\\s*\"+([^\"]+)')[, 2],
    dataset_id = paste0("dataset_name-", extracted_name, "_seed-", extracted_seed)
  ) %>%
  filter(!is.na(extracted_name)) %>%
  mutate(model = recode(module, !!!model_map)) %>%
  group_by(dataset_id, model) %>%
  summarise(mean_time_sec = mean(s, na.rm = TRUE), .groups = "drop")

df_general_meta <- df_metadata %>%
  left_join(marker_map_df, by = "dataset_id") %>%
  mutate(markers_per_pop = n_markers / n_populations) %>%
  filter(!str_detect(dataset_id, regex("sub-sampling", ignore_case = TRUE))) %>%
  filter(!str_detect(dataset_id, regex("Levine", ignore_case = TRUE))) %>%
  mutate(dataset_clean = recode(dataset_id, !!!name_map)) %>%
  filter(!is.na(dataset_clean))

df_plot_f1 <- df_f1_avg %>%
  inner_join(df_general_meta, by = c("dataset" = "dataset_id")) %>%
  filter(!str_detect(model, "Random"))

df_plot_time <- df_time_avg %>%
  inner_join(df_general_meta, by = "dataset_id") %>%
  filter(!str_detect(model, "Random"))

extract_size <- function(name) as.numeric(str_extract(name, "(?<=sub-sampling-)\\d+"))

df_levine_time <- df_perf_raw %>%
  filter(str_detect(params, "sub-sampling")) %>%
  mutate(
    train_size = as.numeric(str_match(params, 'sub-sampling[^:]+:\\s*\"+([^\"]+)')[, 2]),
    model = recode(module, !!!model_map)
  ) %>%
  filter(model %in% names(tool_colors)) %>%
  group_by(model, train_size) %>%
  summarise(mean_time = mean(s, na.rm = TRUE), .groups = "drop")

df_levine_f1 <- df_f1_raw %>%
  filter(str_detect(dataset, "sub-sampling")) %>%
  mutate(
    train_size = extract_size(dataset),
    model = recode(model, !!!model_map)
  ) %>%
  filter(model %in% names(tool_colors)) %>%
  group_by(model, train_size) %>%
  summarise(mean_f1 = mean(f1_weighted_median, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------------------------
# 3. FORCE IDENTICAL LEVELS (prevents weird legend behavior)
# ------------------------------------------------------------------------------
model_levels <- setdiff(names(tool_colors), "Random")
df_plot_f1     <- df_plot_f1     %>% mutate(model = factor(model, levels = model_levels))
df_plot_time   <- df_plot_time   %>% mutate(model = factor(model, levels = model_levels))
df_levine_time <- df_levine_time %>% mutate(model = factor(model, levels = model_levels))
df_levine_f1   <- df_levine_f1   %>% mutate(model = factor(model, levels = model_levels))

my_color_scale <- scale_colour_manual(values = tool_colors[model_levels], name = "Method", drop = FALSE)

legend_theme <- theme(
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.title = element_text(face = "bold", size = 10),
  legend.text = element_text(size = 9),
  legend.key.size = unit(5, "mm"),
  legend.spacing.x = unit(2, "mm")
)

# ------------------------------------------------------------------------------
# 4. GENERATE PLOTS (NO fill MAPPING -> ONLY ONE LEGEND)
# ------------------------------------------------------------------------------

base_point <- geom_point(shape = 16, size = 2.2, alpha = 0.9)
base_line  <- geom_line(linewidth = 0.5, alpha = 0.5)

p_a <- ggplot(df_plot_f1, aes(x = mean_cells, y = mean_f1_weighted, colour = model, group = model)) +
  base_line + base_point +
  my_color_scale +
  scale_x_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Mean Cells / Sample", y = "Mean F1-Score") +
  theme_gb_scatter

p_b <- ggplot(df_plot_f1, aes(x = markers_per_pop, y = mean_f1_weighted, colour = model, group = model)) +
  base_line + base_point +
  my_color_scale +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Number of Markers / Population", y = "Mean F1-Score") +
  theme_gb_scatter

df_plot_time_cells <- df_plot_time %>% arrange(mean_cells)
p_c <- ggplot(df_plot_time_cells, aes(x = mean_cells, y = mean_time_sec, colour = model, group = model)) +
  base_line + base_point +
  my_color_scale +
  scale_x_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
  scale_y_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
  labs(x = "Mean Cells / Sample", y = "Mean Runtime (s)") +
  theme_gb_scatter

df_plot_time_markers <- df_plot_time %>% arrange(markers_per_pop)
p_d <- ggplot(df_plot_time_markers, aes(x = markers_per_pop, y = mean_time_sec, colour = model, group = model)) +
  base_line + base_point +
  my_color_scale +
  scale_y_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
  labs(x = "Number of Markers / Population", y = "Mean Runtime (s)") +
  theme_gb_scatter

p_e <- ggplot(df_levine_time, aes(x = train_size, y = mean_time, colour = model, group = model)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  base_point +
  my_color_scale +
  scale_y_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
  scale_x_continuous(labels = label_number(suffix = "K", scale = 1e-3), breaks = c(10000, 50000, 100000)) +
  labs(x = "Training Set Size", y = "Runtime (s)") +
  theme_gb_scatter

p_f <- ggplot(df_levine_f1, aes(x = train_size, y = mean_f1, colour = model, group = model)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  base_point +
  my_color_scale +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(labels = label_number(suffix = "K", scale = 1e-3), breaks = c(10000, 50000, 100000)) +
  labs(x = "Training Set Size", y = "Mean F1 (Weighted)") +
  theme_gb_scatter

# ------------------------------------------------------------------------------
# 5. BUILD GRID WITH COWPLOT (one shared legend)
# ------------------------------------------------------------------------------

# Extract legend from ONE plot (all have same scale/levels)
legend <- get_legend(p_a + legend_theme)

# Remove legends from all panels
p_a_nl <- p_a + theme(legend.position = "none")
p_b_nl <- p_b + theme(legend.position = "none")
p_c_nl <- p_c + theme(legend.position = "none")
p_d_nl <- p_d + theme(legend.position = "none")
p_e_nl <- p_e + theme(legend.position = "none")
p_f_nl <- p_f + theme(legend.position = "none")

# Add panel tags a-f
p_a_nl <- p_a_nl + labs(tag = "a")
p_b_nl <- p_b_nl + labs(tag = "b")
p_c_nl <- p_c_nl + labs(tag = "c")
p_d_nl <- p_d_nl + labs(tag = "d")
p_e_nl <- p_e_nl + labs(tag = "e")
p_f_nl <- p_f_nl + labs(tag = "f")

tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 12),
  plot.tag.position = c(0, 1),
  plot.tag.justification = c(0, 1)
)

p_a_nl <- p_a_nl + tag_theme
p_b_nl <- p_b_nl + tag_theme
p_c_nl <- p_c_nl + tag_theme
p_d_nl <- p_d_nl + tag_theme
p_e_nl <- p_e_nl + tag_theme
p_f_nl <- p_f_nl + tag_theme

grid <- plot_grid(
  p_a_nl, p_b_nl,
  p_c_nl, p_d_nl,
  p_e_nl, p_f_nl,
  ncol = 2,
  align = "hv"
)

final_fig <- plot_grid(
  grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# ------------------------------------------------------------------------------
# 6. SAVE
# ------------------------------------------------------------------------------
if (!str_detect(output_file, "\\.png$")) {
  output_file <- paste0(output_file, ".png")
}

ggsave(
  output_file,
  plot = final_fig,
  width = 180,
  height = 240,
  units = "mm",
  dpi = 600
)

print(paste("Saved Combined Figure to", output_file))

### Usage example:
# Rscript combined_plot_3.r \
#   --f1_weighted ../ob-blob-metrics/out/metric_collectors/metrics_report/samples_vs_f1_weighted.tsv \
#   --perf_input ../ob-blob-metrics/out/performances.tsv \
#   --meta_json ../ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
#   --output ../ob-pipeline-plots/Figure3_Combined_Grid.png