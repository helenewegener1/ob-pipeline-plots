#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(stringr)
  library(viridis)
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--f1_input"), type = "character", default = NULL, 
              help = "Path to the input f1_macro TSV file",metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "Figure1_Heatmap-boxplot.png", 
              help = "Path to the output PNG file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$f1_input)) {
  print_help(opt_parser)
  stop("Input file (-i / --f1_input) must be provided.", call. = FALSE)
}

input_file <- opt$f1_input
output_file <- opt$output

# ------------------------------------------------------------------------------
# 2. DATA PROCESSING & MAPPINGS
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
  "dataset_name-FR-FCM-ZZRQ_seed-42" = "DCI"
)

model_map <- c(
'cyanno' = "CyAnno",
'cygate' = "CyGATE",
'dgcytof' = "DGCytof",
'gatemeclass' = "GateMeClass",
'lda' = "CyTOF LC",
'random' = "Random"
)

tissue_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 'PBMCs',
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 'PBMCs',
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 'PBMCs',
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 'Stimulated Blood',
  'dataset_name-Samusik_seed-42'                    = 'Bone Marrow',
  'dataset_name-Transformed_seed-42'                = 'PBMCs',
  'dataset_name-flowcyt_seed-42' = 'Bone Marrow',
  'dataset_name-Levine_seed-42' = 'Bone Marrow',
  "dataset_name-FR-FCM-ZZRQ_seed-42" = "PBMCs"
)

populations_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 16,
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 6,
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 6,
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 28,
  'dataset_name-Samusik_seed-42'                    = 15,
  'dataset_name-Transformed_seed-42'                = 8,
  'dataset_name-flowcyt_seed-42' = 5,
  'dataset_name-Levine_seed-42' = 15,
  "dataset_name-FR-FCM-ZZRQ_seed-42" = 8
)

marker_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 37,
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 24,
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 24,
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 38,
  'dataset_name-Samusik_seed-42'                    = 39,
  'dataset_name-Transformed_seed-42'                = 33,
  'dataset_name-flowcyt_seed-42' = 12,
  'dataset_name-Levine_seed-42' = 32,
  "dataset_name-FR-FCM-ZZRQ_seed-42" = 9
)

platform_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = "CyTOF",
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = "FlowCyt",
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = "FlowCyt",
  'dataset_name-FR-FCM-Z3YR_seed-42'                = "CyTOF",
  'dataset_name-Samusik_seed-42'                    = "CyTOF",
  'dataset_name-Transformed_seed-42'                = "CyTOF",
  'dataset_name-flowcyt_seed-42' = "FlowCyt",
  'dataset_name-Levine_seed-42' = 'CyTOF',
  "dataset_name-FR-FCM-ZZRQ_seed-42" = "FlowCyt"
)

df <- read_tsv(input_file, show_col_types = FALSE)

df_clean <- df %>%
  # EXCLUDE subsampled DATASETS
  filter(!str_detect(dataset, regex("sub-sampling", ignore_case = TRUE))) %>%
  mutate(
    platform = platform_map[dataset],
    markers  = marker_map[dataset],
    tissue = tissue_map[dataset],
    model = model_map[model],
    populations = populations_map[dataset],
    clean_name = recode(dataset, !!!name_map),
    display_name = ifelse(!is.na(markers), 
                          paste0(clean_name, " (",tissue,"- M:", markers, "- P:",populations, ")"), 
                          clean_name),
    platform = ifelse(is.na(platform), "Other", platform)
  )

# Filter for Bottom Plot (no Random)
df_no_random <- df_clean %>%
  filter(!str_detect(model, regex("random", ignore_case = TRUE)))

# ------------------------------------------------------------------------------
# 3. DYNAMIC ORDERING (FACET-AWARE)
# ------------------------------------------------------------------------------

dataset_order_df <- df_no_random %>%
  group_by(platform, display_name) %>%
  summarise(global_mean = mean(f1_macro, na.rm = TRUE), .groups = "drop") %>%
  arrange(platform, global_mean)

dataset_order <- dataset_order_df$display_name

tool_order <- df_clean %>%
  group_by(model) %>%
  summarise(global_mean = mean(f1_macro, na.rm = TRUE)) %>%
  arrange(global_mean) %>%
  pull(model)

df_clean$model <- factor(df_clean$model, levels = tool_order)
df_clean$display_name <- factor(df_clean$display_name, levels = dataset_order)
df_no_random$display_name <- factor(df_no_random$display_name, levels = dataset_order)

heatmap_data <- df_clean %>%
  group_by(model, display_name, platform) %>%
  summarise(mean_f1 = mean(f1_macro, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------------------------
# 4. PLOTTING
# ------------------------------------------------------------------------------

theme_gb_base <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    plot.margin = margin(2, 2, 2, 2)
  )

# --- A. HEATMAP WITH PLATFORM FACETING ---
p_heatmap <- ggplot(heatmap_data, aes(x = display_name, y = model, fill = mean_f1)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", mean_f1)), size = 1.8, color = "black") +
  
  facet_grid(. ~ platform, scales = "free_x", space = "free_x") +
  
  scale_fill_viridis_c(option = "viridis", limits = c(0, 1), guide = "none") +
  
  # Removed str_wrap to keep names linear and clean
  scale_x_discrete(position = "top", expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  
  labs(x = "", y = NULL) +
  
  theme_gb_base +
  theme(
    # Angle at 45 or 90 for better readability of long names
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    # Add space between the labels and the "Training set" title
    axis.title.x = element_text(margin = margin(b = 20), size = 9), 
    panel.spacing = unit(8, "pt"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid = element_blank()
  )

# --- B. RIGHT BOXPLOT (No Red Dot) ---
p_right <- ggplot(df_clean, aes(x = f1_macro, y = model)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.25, fill = "white") +
  # stat_summary REMOVED here
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "F1-score", y = NULL) +
  theme_gb_base +
  theme(
    panel.background = element_rect(fill = "grey95", color = NA),
    axis.text.y = element_blank(),
    plot.margin = margin(l = 2)
  )

# --- C. BOTTOM BOXPLOT (No Red Dot, Faceted) ---
p_bottom <- ggplot(df_no_random, aes(x = display_name, y = f1_macro)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.25, fill = "white") +
  # stat_summary REMOVED here
  facet_grid(. ~ platform, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "F1-score") +
  theme_gb_base +
  theme(
    panel.background = element_rect(fill = "grey95", color = NA),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(6, "pt"),
    plot.margin = margin(t = 2)
  )

# ------------------------------------------------------------------------------
# 5. COMPOSITION & SAVING
# ------------------------------------------------------------------------------

design <- "
AAAAB
AAAAB
AAAAB
CCCC#
"

final_plot <- wrap_plots(A = p_heatmap, B = p_right, C = p_bottom, design = design)

# Dynamic Sizing (mm)
n_tools <- length(unique(df_clean$model))
n_datasets <- length(unique(df_clean$display_name))

final_h <- min(max((n_tools * 8) + 80, 120), 225)
final_w <- min(max((n_datasets * 15) + 80, 140), 180)

ggsave(output_file, plot = final_plot, width = final_w, height = final_h, units = "mm", dpi = 600)

print(paste("Saved to", output_file))

###
# Usage example:
# Rscript Fig1_plot.r \
#  --f1_input ./ob-blob-metrics/out/metric_collectors/metrics_report/f1_macro_by_crossvalidation.tsv \
#  --output ./ob-pipeline-plots/Figure1_Heatmap-boxplot.png
###