#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 1. LIBRARIES
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
  library(cowplot)
})

# ------------------------------------------------------------------------------
# 2. COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-f", "--f1_macro"), type = "character", default = NULL,
              help = "Path to f1_macro_by_crossvalidation.tsv", metavar = "character"),
  make_option(c("-p", "--perf_input"), type = "character", default = NULL,
              help = "Path to run_metrics.tsv", metavar = "character"),
  make_option(c("-j", "--meta_json"), type = "character", default = NULL,
              help = "Path to dataset_metadata.json", metavar = "character"),
  # CHANGED: Replaced single output with two distinct output options
  make_option(c("-o", "--out_fig1"), type = "character", default = "Figure_Markers_Pops.png",
              help = "Output filename for Markers and Populations grid", metavar = "character"),
  make_option(c("-q", "--out_fig2"), type = "character", default = "Figure_Cells_TrainSize.png",
              help = "Output filename for Mean Cells and Train Size grid", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


f1_macro <-  "Documents/courses/Benchmarking/repos/ob-blob-metrics/out/metric_collectors/metrics_report/f1_macro_by_crossvalidation.tsv" 
perf_input <- "Documents/courses/Benchmarking/repos/ob-blob-metrics/out/metric_collectors/metrics_report/run_metrics.tsv"
meta_json <- "Documents/courses/Benchmarking/repos/ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json"

out_fig1 <- "Documents/courses/Benchmarking/repos/ob-pipeline-plots/NewNames_Figure3a_Markers_Pops.png"
out_fig2 <- "Documents/courses/Benchmarking/repos/ob-pipeline-plots/NewNames_Cells_TrainSize.png"


# ------------------------------------------------------------------------------
# 3. GLOBAL MAPPINGS
# ------------------------------------------------------------------------------
# name_map <- c(
#   'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 'CV',
#   'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 'HT',
#   'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 'COVID',
#   'dataset_name-FR-FCM-Z3YR_seed-42'                = 'SB',
#   'dataset_name-Samusik_seed-42'                    = 'MBM',
#   'dataset_name-Transformed_seed-42'                = 'TF',
#   'dataset_name-flowcyt_seed-42'                    = 'HBM',
#   'dataset_name-Levine_seed-42'                     = 'LV'
# )

# name_map <- c(
#   'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 'ChikVirusPBMC_Cyt',
#   'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 'PBMC_flow',
#   'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 'covidPBMC_flow',
#   'dataset_name-FR-FCM-Z3YR_seed-42'                = 'StimBlood_Cyt',
#   'dataset_name-Samusik_seed-42'                    = 'MouseBoneMarrow_Cyt',
#   'dataset_name-Transformed_seed-42'                = 'PBMC_Cyt',
#   'dataset_name-flowcyt_seed-42'                    = 'HumanBoneMarrow_flow',
#   'dataset_name-Levine_seed-42'                     = 'HumanBoneMarrow_Cyt',
#   "dataset_name-panel_CD20_seed-42"                 = "DCI-CD20",
#   "dataset_name-panel_CD56_seed-42"                 = "DCI-CD56"
# )

name_map <- c(
  'dataset_name-FR-FCM-Z238_seed-42' = 'ChikVirusPBMC_Cyt',
  'dataset_name-FR-FCM-Z2KP-healthy_seed-42'  = 'PBMC_flow',
  'dataset_name-FR-FCM-Z2KP-covid_seed-42'    = 'covidPBMC_flow',
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 'StimBlood_Cyt',
  'dataset_name-Samusik_seed-42'                    = 'MouseBoneMarrow_Cyt',
  'dataset_name-BodenmillerXL_seed-42'                = 'PBMC_Cyt',
  'dataset_name-FlowCyt_seed-42'                    = 'HumanBoneMarrow_flow',
  'dataset_name-Levine_seed-42'                     = 'HumanBoneMarrow_Cyt'
  # "dataset_name-panel_CD20_seed-42"                 = "DCI-CD20",
  # "dataset_name-panel_CD56_seed-42"                 = "DCI-CD56"
)

populations_map <- c(
  'dataset_name-FR-FCM-Z238_seed-42'  = 16,
  'dataset_name-FR-FCM-Z2KP-healthy_seed-42'  = 6,
  'dataset_name-FR-FCM-Z2KP-covid_seed-42'    = 6,
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 28,
  'dataset_name-Samusik_seed-42'                    = 15,
  'dataset_name-BodenmillerXL_seed-42'                = 8,
  'dataset_name-FlowCyt_seed-42'                    = 5,
  'dataset_name-Levine_seed-42'                     = 15
)

marker_map <- c(
  'dataset_name-FR-FCM-Z238_seed-42' = 37,
  'dataset_name-FR-FCM-Z2KP-healthy_seed-42'  = 24,
  'dataset_name-FR-FCM-Z2KP-covid_seed-42'    = 24,
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 38,
  'dataset_name-Samusik_seed-42'                    = 39,
  'dataset_name-BodenmillerXL_seed-42' = 33,
  'dataset_name-FlowCyt_seed-42' = 12,
  'dataset_name-Levine_seed-42' = 32
)

model_map <- c(
'cyanno' = "CyAnno",
'cygate' = "CyGATE",
'dgcytof' = "DGCyTOF",
'gatemeclass[V]' = "GateMeClass-V",
'gatemeclass[E]' = "GateMeClass-E",
'lda' = "CyTOF LC",
'knn' = "kNN"
)

# --- DEFINING BOLD, "SERIOUS" COLORS (Set1 Palette - No Pastels) ---
tool_colors <- c(
  "CyAnno"        = "#E41A1C",  
  "CyGATE"        = "#377EB8",  
  "DGCyTOF"       = "#4DAF4A",  
  "CyTOF LC"      = "#FF7F00",  
  "GateMeClass-E" = "#6A3D9A",  
  "GateMeClass-V" = "#CAB2D6",  
  "kNN"           = "#17BECF"
)

theme_gb_scatter <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom")

# ------------------------------------------------------------------------------
# 4. DATA PROCESSING
# ------------------------------------------------------------------------------
message("Loading data...")
# df_f1_raw     <- read_tsv(opt$f1_macro, show_col_types = FALSE)
# df_perf_raw   <- read_tsv(opt$perf_input, show_col_types = FALSE)
# metadata_json <- fromJSON(opt$meta_json)

df_f1_raw     <- read_tsv(f1_macro, show_col_types = FALSE)
df_perf_raw   <- read_tsv(perf_input, show_col_types = FALSE)
metadata_json <- fromJSON(meta_json)

# Metadata Setup
cell_stats <- do.call(rbind, lapply(names(metadata_json), function(id) {
  cells <- as.numeric(unlist(metadata_json[[id]]$cells_per_sample))
  data.frame(dataset_id = id, mean_cells = ifelse(all(is.na(cells)), NA, mean(cells, na.rm = TRUE)))
}))

df_general_meta <- data.frame(dataset_id = names(name_map)) %>%
  mutate(n_populations = recode(dataset_id, !!!populations_map),
         n_markers     = recode(dataset_id, !!!marker_map)) %>%
  left_join(cell_stats, by = "dataset_id")

# --- F1 Score Processing ---
df_f1_clean <- df_f1_raw %>%
  filter(!is.na(f1_macro), f1_macro > 0) %>%
  filter(!str_detect(dataset, "panel_CD20|panel_CD56")) %>%
  filter(model != "random")

df_f1_avg <- df_f1_clean %>%
  mutate(model = recode(model, !!!model_map)) %>%
  group_by(dataset, model) %>%
  summarise(mean_f1_macro = mean(f1_macro, na.rm = TRUE), .groups = "drop")

# --- Performance Processing ---
df_perf_processed <- df_perf_raw %>%
  rename(dataset_id = dataset, model_raw = model) %>%
  filter(model_raw != "random")

# Synchronize Performance with F1
df_time_avg <- df_perf_processed %>%
  semi_join(df_f1_clean, by = c("dataset_id" = "dataset", "model_raw" = "model")) %>%
  mutate(model = recode(model_raw, !!!model_map)) %>%
  group_by(dataset_id, model) %>%
  summarise(mean_time_sec = mean(runtime_seconds, na.rm = TRUE), .groups = "drop")

# Join for plotting
df_plot_f1 <- df_f1_avg %>% inner_join(df_general_meta, by = c("dataset" = "dataset_id"))
df_plot_time <- df_time_avg %>% inner_join(df_general_meta, by = "dataset_id")

# --- Sub-sampling Data ---
message("Processing sub-sampling data...")
df_subsampling_f1 <- df_f1_clean %>%
  filter(str_detect(dataset, "sub-sampling")) %>%
  mutate(train_size = as.numeric(str_extract(dataset, "(?<=sub-sampling-)\\d+")),
         model = recode(model, !!!model_map)) %>%
  group_by(model, train_size) %>% 
  summarise(mean_f1 = mean(f1_macro, na.rm = TRUE), .groups = "drop")

df_subsampling_time <- df_perf_processed %>%
  filter(str_detect(dataset_id, "sub-sampling")) %>%
  mutate(train_size = as.numeric(str_extract(dataset_id, "(?<=sub-sampling-)\\d+")),
         model = recode(model_raw, !!!model_map)) %>%
  group_by(model, train_size) %>% 
  summarise(mean_time = mean(runtime_seconds, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------------------------
# 5. PLOTTING
# ------------------------------------------------------------------------------
my_colors <- scale_color_manual(values = tool_colors, name = "Method")
y_f1_limit <- scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

# --- Group 1: Markers and Populations ---
p_markers_f1 <- ggplot(df_plot_f1, aes(n_markers, mean_f1_macro, color=model)) + 
  geom_line(alpha=0.3) + geom_point() + 
  my_colors + y_f1_limit + labs(x="Number of Markers", y="Mean F1") + theme_gb_scatter

p_markers_time <- ggplot(df_plot_time, aes(n_markers, mean_time_sec, color=model)) + 
  geom_line(alpha=0.3) + geom_point() + scale_y_log10() +
  my_colors + labs(x="Number of Markers", y="Mean Runtime (s)") + theme_gb_scatter

p_pops_f1 <- ggplot(df_plot_f1, aes(n_populations, mean_f1_macro, color=model)) + 
  geom_line(alpha=0.3) + geom_point() + 
  my_colors + y_f1_limit + labs(x="Number of Populations", y="Mean F1") + theme_gb_scatter

p_pops_time <- ggplot(df_plot_time, aes(n_populations, mean_time_sec, color=model)) + 
  geom_line(alpha=0.3) + geom_point() + scale_y_log10() +
  my_colors + labs(x="Number of Populations", y="Mean Runtime (s)") + theme_gb_scatter

# --- Group 2: Cells and Train Size ---
p_cells_f1 <- ggplot(df_plot_f1, aes(mean_cells, mean_f1_macro, color=model)) + 
  geom_line(alpha=0.3) + geom_point() + scale_x_log10(labels = label_log()) + 
  my_colors + y_f1_limit + labs(x="Mean Cells per Sample", y="Mean F1") + theme_gb_scatter

p_cells_time <- ggplot(df_plot_time, aes(mean_cells, mean_time_sec, color=model)) + 
  geom_line(alpha=0.3) + geom_point() + scale_x_log10(labels = label_log()) + scale_y_log10() +
  my_colors + labs(x="Mean Cells per Sample", y="Mean Runtime (s)") + theme_gb_scatter

p_train_f1 <- ggplot(df_subsampling_f1, aes(train_size, mean_f1, color=model)) + 
  geom_line() + geom_point() + scale_x_continuous(labels = label_number(suffix="K", scale=1e-3)) +
  my_colors + y_f1_limit + labs(x="Train Size", y="Mean F1") + theme_gb_scatter

p_train_time <- ggplot(df_subsampling_time, aes(train_size, mean_time, color=model)) + 
  geom_line() + geom_point() + scale_y_log10() + scale_x_continuous(labels = label_number(suffix="K", scale=1e-3)) +
  my_colors + labs(x="Train Size", y="Mean Runtime (s)") + theme_gb_scatter

# ------------------------------------------------------------------------------
# 6. GRID & SAVE
# ------------------------------------------------------------------------------
legend <- get_legend(p_cells_f1)

# Helper function to remove legend from individual plots before gridding
strip_legend <- function(p) p + theme(legend.position="none")

# --- GRID 1: Markers & Populations ---
grid1 <- plot_grid(
  strip_legend(p_markers_f1), strip_legend(p_markers_time),
  strip_legend(p_pops_f1),    strip_legend(p_pops_time),
  ncol = 2, labels = "auto", align = "hv"
)
final1 <- plot_grid(grid1, legend, ncol = 1, rel_heights = c(1, 0.1))

# --- GRID 2: Cells & Train Size ---
grid2 <- plot_grid(
  strip_legend(p_cells_f1), strip_legend(p_cells_time),
  strip_legend(p_train_f1), strip_legend(p_train_time),
  ncol = 2, labels = "auto", align = "hv"
)
final2 <- plot_grid(grid2, legend, ncol = 1, rel_heights = c(1, 0.1))

# Adjusted height to 160mm since these are now 2-row grids instead of 4-row
# ggsave(opt$out_fig1, final1, width = 180, height = 160, units = "mm", dpi = 600)
# ggsave(opt$out_fig2, final2, width = 180, height = 160, units = "mm", dpi = 600)

ggsave(out_fig1, final1, width = 180, height = 160, units = "mm", dpi = 600)
ggsave(out_fig2, final2, width = 180, height = 160, units = "mm", dpi = 600)


message(paste("Figure 1 saved to:", opt$out_fig1))
message(paste("Figure 2 saved to:", opt$out_fig2))

### Usage example:
#Rscript ./combined_plot_3.r \
#  --f1_macro ../ob-blob-metrics/out/metric_collectors/metrics_report/f1_macro_by_crossvalidation.tsv \
#  --perf_input ../ob-blob-metrics/out/metric_collectors/metrics_report/run_metrics.tsv \
#  --meta_json ../ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
#  --out_fig1 ../ob-pipeline-plots/Figure3a_Markers_Pops.png \
#  --out_fig2 ../ob-pipeline-plots/Figure3b_Cells_traSize.png