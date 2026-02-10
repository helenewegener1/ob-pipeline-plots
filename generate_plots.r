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
  library(scales)
  library(jsonlite)
  library(grid)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--f1_cv"), type = "character", help = "f1_macro_by_crossvalidation.tsv"),
  make_option(c("--f1_weighted"), type = "character", help = "samples_vs_f1_weighted.tsv"),
  make_option(c("--conf_matrix"), type = "character", help = "per_population_confusion.tsv"),
  make_option(c("--perf_tsv"), type = "character", help = "performances.tsv"),
  make_option(c("--meta_json"), type = "character", help = "dataset_metadata.json"),
  make_option(c("-o", "--out_dir"), type = "character", default = "./plots", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Ensure output directory exists
if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# SHARED MAPPINGS & THEMES
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

tool_colors <- c(
  "CyAnno"      = "#E41A1C",  # Red
  "CyGATE"      = "#377EB8",  # Blue
  "DGCytof"     = "#4DAF4A",  # Green
  "GateMeClass" = "#984EA3",  # Purple
  "CyTOF LC"    = "#FF7F00",  # Orange
  "Random"      = "#525252"   # Dark Grey
)

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

# ------------------------------------------------------------------------------
# FIGURE 1: HEATMAP & COMPOSITE BOXPLOTS
# ------------------------------------------------------------------------------
cat("Generating Figure 1 (Composite Heatmap)...\n")
df1 <- read_tsv(opt$f1_cv, show_col_types = FALSE)

df1_clean <- df1 %>%
  # EXCLUDE subsampled DATASETS
  filter(!str_detect(dataset, regex("sub-sampling", ignore_case = TRUE))) %>%
  mutate(
    platform    = platform_map[dataset],
    markers     = marker_map[dataset],
    tissue      = tissue_map[dataset],
    model       = model_map[model],
    populations = populations_map[dataset],
    clean_name  = recode(dataset, !!!name_map),
    # UPDATED DISPLAY NAME LOGIC
    display_name = ifelse(!is.na(markers), 
                          paste0(clean_name, " (",tissue,"- M:", markers, "- P:",populations, ")"), 
                          clean_name),
    platform = ifelse(is.na(platform), "Other", platform)
  )

# Filter for Bottom Plot (no Random)
df1_no_random <- df1_clean %>%
  filter(!str_detect(model, regex("random", ignore_case = TRUE)))

# --- DYNAMIC ORDERING ---

dataset_order_df1 <- df1_no_random %>%
  group_by(platform, display_name) %>%
  summarise(global_mean = mean(f1_macro, na.rm = TRUE), .groups = "drop") %>%
  arrange(platform, global_mean)

dataset1_order <- dataset_order_df1$display_name

tool_order1 <- df1_clean %>%
  group_by(model) %>%
  summarise(global_mean = mean(f1_macro, na.rm = TRUE), .groups = "drop") %>%
  arrange(global_mean) %>%
  pull(model)

df1_clean$model <- factor(df1_clean$model, levels = tool_order1)
df1_clean$display_name <- factor(df1_clean$display_name, levels = dataset1_order)
df1_no_random$display_name <- factor(df1_no_random$display_name, levels = dataset1_order)

heatmap_data1 <- df1_clean %>%
  group_by(model, display_name, platform) %>%
  summarise(mean_f1 = mean(f1_macro, na.rm = TRUE), .groups = "drop")

# --- PLOTTING ---

theme_gb_base1 <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    plot.margin = margin(2, 2, 2, 2)
  )

# A. HEATMAP
p1_heatmap <- ggplot(heatmap_data1, aes(x = display_name, y = model, fill = mean_f1)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", mean_f1)), size = 1.8, color = "black") +
  facet_grid(. ~ platform, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(option = "viridis", limits = c(0, 1), guide = "none") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "", y = NULL) +
  theme_gb_base1 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.title.x = element_text(margin = margin(b = 20), size = 9), 
    panel.spacing = unit(8, "pt"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid = element_blank()
  )

# B. RIGHT BOXPLOT
p1_right <- ggplot(df1_clean, aes(x = f1_macro, y = model)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.25, fill = "white") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "F1-score", y = NULL) +
  theme_gb_base1 +
  theme(
    panel.background = element_rect(fill = "grey95", color = NA),
    axis.text.y = element_blank(),
    plot.margin = margin(l = 2)
  )

# C. BOTTOM BOXPLOT
p1_bottom <- ggplot(df1_no_random, aes(x = display_name, y = f1_macro)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.25, fill = "white") +
  facet_grid(. ~ platform, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "F1-score") +
  theme_gb_base1 +
  theme(
    panel.background = element_rect(fill = "grey95", color = NA),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(6, "pt"),
    plot.margin = margin(t = 2)
  )

# --- COMPOSITION & SAVING ---

design1 <- "
AAAAB
AAAAB
AAAAB
CCCC#
"

final_plot1 <- wrap_plots(A = p1_heatmap, B = p1_right, C = p1_bottom, design = design1)

n_tools1 <- length(unique(df1_clean$model))
n_datasets1 <- length(unique(df1_clean$display_name))

final_h1 <- min(max((n_tools1 * 8) + 80, 120), 225)
final_w1 <- min(max((n_datasets1 * 15) + 80, 140), 180)

output_fig1 <- file.path(opt$out_dir, "Figure1_Heatmap_Composite.png")
ggsave(output_fig1, plot = final_plot1, width = final_w1, height = final_h1, units = "mm", dpi = 600)

cat(paste0("Success: Figure 1 saved to ", output_fig1, "\n"))

# ------------------------------------------------------------------------------
# FIGURE 2: PER-POPULATION BOXPLOTS (DISTRIBUTION ACROSS CV RUNS)
# ------------------------------------------------------------------------------
cat("Generating Figure 2 (Individual Dataset Plots)...\n")
df2 <- read_tsv(opt$conf_matrix, show_col_types = FALSE)

# 3. DATA PROCESSING
df2_processed <- df2 %>%
  filter(!grepl("sub-sampling", dataset, ignore.case = TRUE)) %>%
  mutate(
    dataset = recode(dataset, !!!name_map),
    model   = recode(model, !!!model_map) 
  ) %>%
  mutate(
    across(c(tp, fp, fn), as.numeric),
    population_name = basename(trimws(population_name)), 
    f1_score = ifelse((2*tp + fp + fn) > 0, (2*tp) / (2*tp + fp + fn), 0),
    actual_count = tp + fn
  )

# 4. METADATA (Calculations for Percentage and Ranks)
pop_meta2 <- df2_processed %>%
  group_by(dataset, population_name) %>%
  summarize(total_count = first(actual_count), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(
    pct = (total_count / sum(total_count)) * 100,
    abundance_class = ifelse(pct < 1, "Rare (<1%)", "Common (>=1%)"),
    pop_label = paste0(population_name, " (", round(pct, 2), "%)")
  ) %>%
  group_by(dataset, abundance_class) %>%
  mutate(rank_within_class = row_number(desc(pct))) %>% 
  ungroup()

df2_final <- df2_processed %>%
  inner_join(pop_meta2, by = c("dataset", "population_name"))

# 5. THEMES
# Theme 1: Large Individual Plots
gb_theme_full <- theme_bw(base_size = 9) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    panel.grid.major = element_line(color = "grey92"), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 4, r = 4, b = 4, l = 12, unit = "mm"), # Large bottom margin for labels
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color="black"),
    axis.text.y = element_text(size = 8, color="black"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 9, hjust = 0, margin = margin(b=2))
  )

# Theme 2: Compact Grid Subplots
gb_theme_subplot <- theme_bw(base_size = 8) + 
  theme(
    # Clean Borders
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.2), 
    panel.grid.minor = element_blank(),
    
    # Strip (Header) Styling - Clean White
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 7, color = "black"),
    
    # CRITICAL: Margins to prevent name cutting
    # Top=2mm, Right=2mm, Bottom=35mm (SPACE FOR NAMES), Left=2mm
    plot.margin = margin(2, 2, 2, 16, "mm"),
    
    # Axis Text
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5, color="black"),
    axis.text.y = element_text(size = 6, color="black"),
    axis.title = element_blank(),
    
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 9, hjust = 0, margin = margin(b=3))
  )

# 6. GENERATION LOOP
unique_datasets <- unique(df2_final$dataset)
plot_list_for_grid <- list()

for (ds_name in unique_datasets) {
  
  ds_data <- df2_final %>% filter(dataset == ds_name)
  
  # -------------------------------------------------------
  # A. Fig2_Full_[dataset].png (ALL POPULATIONS)
  # -------------------------------------------------------
  p1_data <- ds_data %>% mutate(pop_label = reorder(pop_label, -actual_count))
  
  p1 <- ggplot(p1_data, aes(x = pop_label, y = f1_score, fill = model, color = model)) +
    geom_boxplot(outlier.size = 0.2, lwd = 0.35, alpha = 0.7, width = 0.75, fatten = 1) +
    scale_fill_manual(values = tool_colors, name = "Method") +
    scale_color_manual(values = tool_colors, name = "Method") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = paste(ds_name, "- Full Population Profile"), x = NULL, y = "F1-Score") +
    gb_theme_full
  
  w1 <- max(150, length(unique(p1_data$pop_label)) * 12)
  ggsave(file.path(opt$out_dir, paste0("Fig2_Full_", ds_name, ".png")), 
         p1, width = w1, height = 120, units = "mm", dpi = 600)
  
  # -------------------------------------------------------
  # B. Fig2_Extremes_[dataset].png (TOP 5 COMMON + RARE)
  # -------------------------------------------------------
  p2_data <- ds_data %>%
    filter(abundance_class == "Rare (<1%)" | (abundance_class == "Common (>=1%)" & rank_within_class <= 5)) %>%
    mutate(
      abundance_class = factor(abundance_class, levels = c("Common (>=1%)", "Rare (<1%)")),
      pop_label = reorder(pop_label, -pct)
    )
  
  if(nrow(p2_data) == 0) next
  
  # Save Individual Extremes Plot
  p2_full <- ggplot(p2_data, aes(x = pop_label, y = f1_score, fill = model, color = model)) +
    geom_boxplot(outlier.size = 0.2, lwd = 0.35, alpha = 0.7, width = 0.75) +
    facet_grid(. ~ abundance_class, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = tool_colors, name = "Method") +
    scale_color_manual(values = tool_colors, name = "Method") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 0.5, 1)) +
    labs(title = paste(ds_name, "- Extremes"), x = NULL, y = "F1-Score") +
    gb_theme_full
  
  ggsave(file.path(opt$out_dir, paste0("Fig2_Extremes_", ds_name, ".png")), 
         p2_full, width = 180, height = 120, units = "mm", dpi = 600)
  
  # -------------------------------------------------------
  # C. PREPARE GRID OBJECT
  # -------------------------------------------------------
  has_common <- any(p2_data$abundance_class == "Common (>=1%)")
  has_rare   <- any(p2_data$abundance_class == "Rare (<1%)")
  
  if (has_common & has_rare) {
    p_grid <- ggplot(p2_data, aes(x = pop_label, y = f1_score, fill = model, color = model)) +
      geom_boxplot(outlier.size = 0.1, lwd = 0.25, alpha = 0.7, width = 0.65, fatten = 1) +
      facet_grid(. ~ abundance_class, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = tool_colors) +
      scale_color_manual(values = tool_colors) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      labs(title = ds_name) +
      gb_theme_subplot
    
    plot_list_for_grid[[ds_name]] <- p_grid
  }
}

# 7. ASSEMBLE A4 GRID & SAVE BOTH VERSIONS
if(length(plot_list_for_grid) > 0) {
  
  plot_list_for_grid <- plot_list_for_grid[sort(names(plot_list_for_grid))]
  
  # Legend
  legend_plot <- ggplot(df2_final, aes(x=model, y=f1_score, fill=model, color=model)) +
    geom_boxplot(alpha=0.7) +
    scale_fill_manual(values = tool_colors, name = "Method") +
    scale_color_manual(values = tool_colors, name = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.text = element_text(size=8, face="bold"))
  
  common_legend <- get_legend(legend_plot)
  
  # Grid
  p_grid_final <- plot_grid(plotlist = plot_list_for_grid, ncol = 2, align = 'hv', axis = 'tb')
  
  final_figure <- plot_grid(
    p_grid_final,
    common_legend,
    ncol = 1,
    rel_heights = c(1, 0.05)
  )
  
  # SAVE 1: Figure2_MultiPanel_A4.png
  outfile_grid_1 <- file.path(opt$out_dir, "Figure2_MultiPanel_A4.png")
  ggsave(outfile_grid_1, plot = final_figure, width = 210, height = 297, units = "mm", dpi = 800)
  
  # SAVE 2: Figure2_A4_Top5_Refined.png
  outfile_grid_2 <- file.path(opt$out_dir, "Figure2_A4_Top5_Refined.png")
  ggsave(outfile_grid_2, plot = final_figure, width = 210, height = 297, units = "mm", dpi = 800)
  
  cat(paste("Success: Saved Full, Extremes, and Grids to:", opt$out_dir, "\n"))
  
} else {
  cat("No datasets qualified for the Grid (need both Common & Rare).\n")
}

# ------------------------------------------------------------------------------
# SHARED METADATA EXTRACTION & THEME FOR FIG 3 & 4
# ------------------------------------------------------------------------------
cat("Processing Metadata for Figures 3 & 4...\n")
metadata_json <- fromJSON(opt$meta_json)

extract_metadata_stats <- function(meta_list) {
  ds_ids <- names(meta_list)
  ds_ids_filtered <- ds_ids[!str_detect(ds_ids, regex("Levine", ignore_case = TRUE))]
  
  results <- lapply(ds_ids_filtered, function(id) {
    entry <- meta_list[[id]]
    cells <- as.numeric(entry$cells_per_sample)
    avg_cells <- if(all(is.na(cells))) NA else mean(cells, na.rm = TRUE)
    
    data.frame(
      dataset_id    = id,
      n_markers     = if(id %in% names(marker_map)) as.numeric(marker_map[id]) else NA,
      n_samples     = as.numeric(entry$sample_count),
      n_populations = as.numeric(entry$population_count),
      mean_cells    = avg_cells,
      stringsAsFactors = FALSE
    )
  })
  return(do.call(rbind, results))
}

df_metadata_shared <- extract_metadata_stats(metadata_json)

theme_gb_scatter_shared <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    panel.background = element_rect(fill = "grey98", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# FIGURE 3: PERFORMANCE vs METADATA (ACCURACY)
# ------------------------------------------------------------------------------
cat("Generating Figure 3 (Metadata vs Performance)...\n")

# 1. LOAD DATA & SETUP
df_f1_raw <- read_tsv(opt$f1_weighted, show_col_types = FALSE)

# Ensure Set1 Colors are defined (if not already global)
tool_colors <- c(
  "CyAnno"      = "#E41A1C",
  "CyGATE"      = "#377EB8",
  "DGCytof"     = "#4DAF4A",
  "GateMeClass" = "#984EA3",
  "CyTOF LC"    = "#FF7F00",  # Orange
  "Random"      = "#525252"
)

# Hardcoded Marker Map (Fixes NAs for specific datasets)
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

# 2. DATA PROCESSING

# A. Summarize F1
df_f1_avg <- df_f1_raw %>%
  mutate(model = recode(model, !!!model_map)) %>% 
  group_by(dataset, model) %>%
  summarise(mean_f1_weighted = mean(f1_weighted_median, na.rm = TRUE), .groups = "drop")

# B. Metadata Extraction (Robust Re-extraction)
# NOTE: This assumes 'metadata_json' object exists in your environment. 
# If 'df_metadata_shared' was used before, we rebuild it here for safety to match your Code A logic.
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

# Generate metadata df (Assuming metadata_json is loaded)
df_metadata <- extract_metadata_stats(metadata_json)

# C. Merge and Final Cleanup
df3_plot <- df_f1_avg %>%
  left_join(df_metadata, by = c("dataset" = "dataset_id")) %>%
  left_join(marker_map_df, by = c("dataset" = "dataset_id")) %>% # Fixes Marker NAs
  filter(!str_detect(dataset, regex("sub-sampling", ignore_case = TRUE))) %>%
  filter(!str_detect(dataset, regex("Levine", ignore_case = TRUE))) %>% 
  filter(!str_detect(model, regex("random", ignore_case = TRUE))) %>%
  mutate(dataset_clean = recode(dataset, !!!name_map)) %>%
  filter(!is.na(dataset_clean))

# 3. PLOTTING FUNCTIONS (Serious Theme)

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

create_fig3_panel <- function(data, x_var, x_label, is_log = FALSE) {
  p <- ggplot(data, aes(x = .data[[x_var]], y = mean_f1_weighted, color = model, fill = model)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, alpha = 0.5) +
    geom_point(shape = 21, color = "black", stroke = 0.2, size = 2.5, alpha = 0.8) +
    
    # Apply Consistent Set1 Colors
    scale_color_manual(values = tool_colors, name = "Method") +
    scale_fill_manual(values = tool_colors, name = "Method") +
    
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = x_label, y = "Mean Weighted F1-Score") +
    theme_gb_scatter
  
  if(is_log) {
    p <- p + scale_x_log10(
      labels = label_log(), 
      breaks = trans_breaks("log10", function(x) 10^x)
    )
  }
  return(p)
}

# 4. GENERATE PANELS
p3_cells   <- create_fig3_panel(df3_plot, "mean_cells", "Mean Cells / Sample", is_log = TRUE)
p3_markers <- create_fig3_panel(df3_plot, "n_markers", "Number of Markers")
p3_samples <- create_fig3_panel(df3_plot, "n_samples", "Number of Samples")
p3_pops    <- create_fig3_panel(df3_plot, "n_populations", "Number of Populations")

# 5. ASSEMBLE
final_fig3 <- (p3_cells + p3_markers) / (p3_samples + p3_pops) + 
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
outfile_fig3 <- file.path(opt$out_dir, "Figure3_Accuracy_Metadata.png")
ggsave(outfile_fig3, final_fig3, width = 180, height = 160, units = "mm", dpi = 600)

cat(paste("Success: Figure 3 saved to", outfile_fig3, "\n"))

# ------------------------------------------------------------------------------
# FIGURE 4: RUNTIME vs METADATA
# ------------------------------------------------------------------------------
cat("Generating Figure 4 (Metadata vs Runtime)...\n")

# 1. LOAD DATA
df4_perf_raw <- read_tsv(opt$perf_tsv, show_col_types = FALSE)

# 2. DATA PROCESSING

# A. Parse Performance Data
df4_time_avg <- df4_perf_raw %>%
  # Exclude internal pipeline steps
  filter(!module %in% c("data_import", "data_preprocessing", "flow_metrics", "f1_score", "metric_collector", "metrics")) %>%
  # Regex to extract dataset name and seed
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

# B. Metadata Extraction (Robust)
# (Re-using the robust function defined in Fig 3 block, or ensuring it runs here)
# If running independently, ensure 'extract_metadata_stats' is defined.
df_metadata_4 <- extract_metadata_stats(metadata_json)

# C. Merge and Final Cleanup
df4_plot <- df4_time_avg %>%
  left_join(df_metadata_4, by = c("dataset_id" = "dataset_id")) %>%
  left_join(marker_map_df, by = c("dataset_id" = "dataset_id")) %>% # Marker Fix
  filter(!str_detect(dataset_id, regex("sub-sampling", ignore_case = TRUE))) %>%
  filter(!str_detect(dataset_id, regex("Levine", ignore_case = TRUE))) %>%
  filter(!str_detect(model, regex("random", ignore_case = TRUE))) %>%
  mutate(dataset_clean = recode(dataset_id, !!!name_map)) %>%
  filter(!is.na(dataset_clean))

# 3. PLOTTING FUNCTIONS

create_fig4_panel <- function(data, x_var, x_label, is_log_x = FALSE) {
  
  p <- ggplot(data, aes(x = .data[[x_var]], y = mean_time_sec, color = model, fill = model)) +
    # Linear Trend Line
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, alpha = 0.5) +
    
    # Scatter Points
    geom_point(shape = 21, color = "black", stroke = 0.2, size = 2.5, alpha = 0.8) +
    
    # Consistent Colors (Set1)
    scale_color_manual(values = tool_colors, name = "Method") +
    scale_fill_manual(values = tool_colors, name = "Method") +
    
    # Y AXIS: LOG SCALE FOR RUNTIME
    scale_y_log10(
      labels = label_log(), 
      breaks = trans_breaks("log10", function(x) 10^x)
    ) +
    
    labs(x = x_label, y = "Mean Runtime (s)") +
    theme_gb_scatter # Re-using the theme defined in Fig 3
  
  if(is_log_x) {
    p <- p + scale_x_log10(
      labels = label_log(), 
      breaks = trans_breaks("log10", function(x) 10^x)
    )
  }
  return(p)
}

# 4. GENERATE PANELS
p4_cells   <- create_fig4_panel(df4_plot, "mean_cells", "Mean Cells / Sample", is_log_x = TRUE)
p4_markers <- create_fig4_panel(df4_plot, "n_markers", "Number of Markers")
p4_samples <- create_fig4_panel(df4_plot, "n_samples", "Number of Samples")
p4_pops    <- create_fig4_panel(df4_plot, "n_populations", "Number of Populations")

# 5. ASSEMBLE
final_fig4 <- (p4_cells + p4_markers) / (p4_samples + p4_pops) + 
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
outfile_fig4 <- file.path(opt$out_dir, "Figure4_Runtime_Metadata.png")
ggsave(outfile_fig4, final_fig4, width = 180, height = 160, units = "mm", dpi = 600)

cat(paste("Success: Figure 4 saved to", outfile_fig4, "\n"))

# ------------------------------------------------------------------------------
# FIGURE 5: SCALABILITY (SUB-SAMPLING)
# ------------------------------------------------------------------------------
cat("Generating Figure 5 (Scalability Panels)...\n")

# 1. Helper for Dataset String Parsing
extract_size <- function(name) {
  as.numeric(str_extract(name, "(?<=sub-sampling-)\\d+"))
}
# 2. DATA PROCESSING

# A. Process Performance (Time)
df5_time <- read_tsv(opt$perf_tsv, show_col_types = FALSE) %>%
  filter(str_detect(params, "sub-sampling")) %>%
  mutate(
    # Extract size from params string
    train_size = as.numeric(str_match(params, 'sub-sampling[^:]+:\\s*\"+([^\"]+)')[,2]),
    # Clean Tool Name
    model = recode(module, !!!model_map)
  ) %>%
  filter(!module %in% c("metrics", "flow_metrics", "data_import", "data_preprocessing", "random")) %>%
  filter(model %in% names(tool_colors)) %>% # Safety filter
  group_by(model, train_size) %>%
  summarise(mean_time = mean(s, na.rm = TRUE), .groups = "drop")

# B. Process F1 (Accuracy)
df5_f1 <- read_tsv(opt$f1_weighted, show_col_types = FALSE) %>%
  filter(str_detect(dataset, "sub-sampling")) %>%
  mutate(
    # Extract size from dataset name
    train_size = extract_size(dataset),
    # Clean Tool Name
    model = recode(model, !!!model_map)
  ) %>%
  filter(!str_detect(model, "Random")) %>%
  filter(model %in% names(tool_colors)) %>%
  group_by(model, train_size) %>%
  summarise(mean_f1 = mean(f1_weighted_median, na.rm = TRUE), .groups = "drop")


# 3. PLOTTING SETUP

theme_gb_line <- theme_bw(base_size = 9) +
  theme(
    text = element_text(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    legend.position = "none", # Collected later
    plot.title = element_text(face = "bold", size = 10)
  )

# 4. GENERATE PLOTS

# Plot A: Time vs Training Size
p5_time <- ggplot(df5_time, aes(x = train_size, y = mean_time, color = model, group = model)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  # Points: Filled with color, thin black border
  geom_point(aes(fill = model), shape = 21, color = "black", stroke = 0.2, size = 2.5) +
  
  # Colors
  scale_color_manual(values = tool_colors, name = "Method") +
  scale_fill_manual(values = tool_colors, name = "Method") +
  
  # Scales
  scale_y_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
  scale_x_continuous(labels = label_number(suffix = "K", scale = 1e-3), breaks = c(10000, 50000, 100000)) +
  
  labs(title = "Computational Scalability", x = NULL, y = "Runtime (s)") +
  theme_gb_line +
  # FORCE 2 ROWS IN LEGEND
  guides(color = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE))

# Plot B: F1 vs Training Size
p5_f1 <- ggplot(df5_f1, aes(x = train_size, y = mean_f1, color = model, group = model)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(aes(fill = model), shape = 21, color = "black", stroke = 0.2, size = 2.5) +
  
  # Colors
  scale_color_manual(values = tool_colors, name = "Method") +
  scale_fill_manual(values = tool_colors, name = "Method") +
  
  # Scales
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(labels = label_number(suffix = "K", scale = 1e-3), breaks = c(10000, 50000, 100000)) +
  
  labs(title = "Classification Stability", x = "Training Set Size", y = "Mean F1 (Weighted)") +
  theme_gb_line +
  # FORCE 2 ROWS IN LEGEND
  guides(color = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE))

# 5. COMBINE & SAVE
final_fig5 <- (p5_time / p5_f1) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a') & 
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 9),
    legend.key.size = unit(5, "mm"),
    plot.tag = element_text(face = "bold", size = 12)
  )

outfile_fig5 <- file.path(opt$out_dir, "Figure5_Healthy_Scalability_Final.png")

ggsave(outfile_fig5, 
       plot = final_fig5, 
       width = 120, 
       height = 150, 
       units = "mm", 
       dpi = 800)

cat(paste("Success: Figure 5 saved to", outfile_fig5, "\n"))
cat(paste0("\n[ALL TASKS COMPLETE] All plots saved to: ", opt$out_dir, "\n"))
###
# Usage example:
#Rscript generate_plots.r \
#  --f1_cv /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/metric_collectors/metrics_report/f1_macro_by_crossvalidation.tsv \
#  --f1_weighted /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/metric_collectors/metrics_report/samples_vs_f1_weighted.tsv \
#  --conf_matrix /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/metric_collectors/metrics_report/per_population_confusion.tsv \
#  --perf_tsv /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/performances.tsv \
#  --meta_json /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
#  --out_dir /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-pipeline-plots
###