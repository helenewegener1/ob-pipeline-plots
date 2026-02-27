#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(viridis)
  library(patchwork)
  library(grid)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--conf_input"), type = "character", default = NULL, 
              help = "Path to the per_population_confusion.tsv file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./", 
              help = "Directory where the individual PNGs will be saved", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$conf_input)) {
  print_help(opt_parser)
  stop("Input file (-i / --conf_input) must be provided.", call. = FALSE)
}

input_file <- opt$conf_input
output_dir <- opt$output_dir

input_file <- "Documents/courses/Benchmarking/repos/ob-blob-metrics/out/metric_collectors/metrics_report/per_population_confusion.tsv"
output_dir <- "Documents/courses/Benchmarking/repos/ob-pipeline-plots/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 1. DATA LOADING & PROCESSING
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
'dgcytof' = "DGCyTOF",
'gatemeclass[V]' = "GateMeClass-V",
'gatemeclass[E]' = "GateMeClass-E",
'lda' = "CyTOF LC",
'knn' = "kNN",
'random' = "Random"
)

# --- DEFINING BOLD, "SERIOUS" COLORS (Set1 Palette - No Pastels) ---
tool_colors <- c(
  "CyAnno"        = "#E41A1C",  # Bold Red
  "CyGATE"        = "#377EB8",  # Strong Blue
  "DGCyTOF"       = "#4DAF4A",  # Vivid Green
  "CyTOF LC"      = "#FF7F00",  # Strong Orange
  "GateMeClass-E" = "#6A3D9A",  # Deep, highly saturated Purple
  "GateMeClass-V" = "#CAB2D6",  # Light Lilac (pairs logically with E, but clearly distinct)
  "kNN"           = "#17BECF",  # Bright Teal/Cyan (moves it entirely away from the purples/reds)
  "Random"        = "#525252"   # Dark Grey
)

df <- read_tsv(input_file, show_col_types = FALSE)

df_processed <- df %>%
  filter(!grepl("sub-sampling", dataset, ignore.case = TRUE)) %>%
  # EXCLUDE DCI DATASETS
  filter(!str_detect(dataset, "panel_CD20|panel_CD56")) %>%
  mutate(
    dataset = recode(dataset, !!!name_map),
    model   = recode(model, !!!model_map) 
  ) %>%
  mutate(
    across(c(tp, fp, fn), as.numeric),
    # Trim whitespace and keep only the actual population name
    population_name = basename(trimws(population_name)), 
    f1_score = ifelse((2*tp + fp + fn) > 0, (2*tp) / (2*tp + fp + fn), 0),
    actual_count = tp + fn
  )

# ------------------------------------------------------------------------------
# 3. METADATA CALCULATION & VALIDATION
# ------------------------------------------------------------------------------
pop_meta <- df_processed %>%
  # 1. Use max() instead of first() to avoid NA/0 errors from failed models
  group_by(dataset, population_name) %>%
  summarize(total_count = max(actual_count, na.rm = TRUE), .groups = "drop") %>%
  
  # 2. Explicitly handle NAs in the sum
  group_by(dataset) %>%
  mutate(
    dataset_total = sum(total_count, na.rm = TRUE),
    pct = (total_count / dataset_total) * 100,
    # abundance_class IF in top 3 "top_3" if bottom 3 "bottom_3" else NA
    abundance_class = case_when(
      min_rank(desc(total_count)) <= 3  ~ "Most prevalent",
      min_rank(total_count)       <= 3  ~ "Least prevalent",
      TRUE                              ~ NA_character_
    ),
    # abundance_class = ifelse(pct < 1, "Rare (<1%)", "Common (>=1%)"),
    pop_label = paste0(population_name, " (", round(pct, 2), "%)"),
    pop_label_newline = paste0(population_name, "\n", round(pct, 2), "%")
  ) %>%
  
  # 3. Use min_rank() which handles ties gracefully (unlike row_number)
  group_by(dataset, abundance_class) %>%
  mutate(rank_within_class = min_rank(desc(total_count))) %>%
  ungroup()

# --- 100% VALIDATION CHECK ---
pct_check <- pop_meta %>%
  group_by(dataset) %>%
  summarize(total_pct = sum(pct), .groups = "drop") %>%
  # We use near() with a small tolerance because floating point math can cause a sum to be 99.9999999
  filter(!dplyr::near(total_pct, 100, tol = 0.01))

if (nrow(pct_check) > 0) {
  # Create a detailed error message listing the problematic datasets
  error_msg <- paste(
    "Error: The following datasets have subpopulation percentages that do not sum to 100%:\n",
    paste(" -", pct_check$dataset, ":", round(pct_check$total_pct, 2), "%", collapse = "\n")
  )
  stop(error_msg, call. = FALSE)
}
# -----------------------------

df_final <- df_processed %>%
  inner_join(pop_meta, by = c("dataset", "population_name")) 

# df_final <- df_processed %>%
#   left_join(pop_meta, by = c("dataset", "population_name")) %>% distinct()

# ------------------------------------------------------------------------------
# 4. THEME & GENERATION LOOP (Full Profile & Extremes)
# ------------------------------------------------------------------------------
gb_theme <- theme_bw(base_size = 9) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey92"), 
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 9, color = "black"),
    
    # Increased bottom margin to 35mm to accommodate long, rotated labels
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm"),
    
    # Text styling
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4, color="black"),
    axis.text.y = element_text(size = 8, color="black"),
    axis.title.y = element_text(face="bold", size=9),
    
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5)
  )

unique_datasets <- unique(df_final$dataset)

for (ds_name in unique_datasets) {
  
  # ds_name <- "CV"
  
  ds_data <- df_final %>% filter(dataset == ds_name) 
  ds_name_n_cells <- df_final %>% filter(dataset == ds_name) %>% select(dataset_total) %>% unique() %>% pull()
  
  pop_meta %>% filter(dataset == ds_name)
  
  # --- PLOT TYPE 1: FULL PROFILE ---
  p1_data <- ds_data %>% mutate(pop_label = reorder(pop_label, -actual_count))
  
  p1 <- ggplot(p1_data, aes(x = pop_label, y = f1_score, fill = model, color = model)) +
    geom_boxplot(
      outlier.size = 0.2, 
      lwd = 0.35,       
      alpha = 0.7,     # Fill transparency
      width = 0.75,
      fatten = 1       # Median line thickness
    ) +
    
    # APPLY BOLD COLORS
    scale_fill_manual(values = tool_colors, name = "Method") +
    scale_color_manual(values = tool_colors, name = "Method") +

    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = paste(ds_name, "- Full Population Profile"), x = NULL, y = "F1-Score", subtitle = paste("N cells:", format(ds_name_n_cells, big.mark = ","))) +
    gb_theme
  
  # Dynamic Width: Ensure at least 15cm, add space for many populations
  w1 <- max(150, length(unique(p1_data$pop_label)) * 12)
  
  ggsave(file.path(output_dir, paste0("Fig2_Full_", ds_name, ".png")), p1, width = w1, height = 120, units = "mm", dpi = 600)

  # --- PLOT TYPE 2: EXTREMES (TOP 5 COMMON VS ALL RARE) ---
  p2_data <- ds_data %>%
    filter(abundance_class == "Most prevalent" | (abundance_class == "Least prevalent" & rank_within_class <= 5)) %>%
    mutate(
      abundance_class = factor(abundance_class, levels = c("Most prevalent", "Least prevalent")),
      pop_label = reorder(pop_label, -actual_count)
    )
  
  n_pop <- p2_data$population_name %>% unique() %>% length()

  p2 <- ggplot(p2_data, aes(x = pop_label, y = f1_score, fill = model, color = model)) +
    geom_boxplot(
      outlier.size = 0.2, 
      lwd = 0.35, 
      alpha = 0.7,
      width = 0.75
    ) +
    
    (if (n_pop >= 6) facet_grid(. ~ abundance_class, scales = "free_x", space = "free_x") else NULL) +
    
    # APPLY BOLD COLORS
    scale_fill_manual(values = tool_colors, name = "Method") +
    scale_color_manual(values = tool_colors, name = "Method") +
    
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = paste(ds_name, "- Most and least prevalent cell populations"), 
      x = NULL, 
      y = "F1-Score",
      subtitle = paste("N cells:", format(ds_name_n_cells, big.mark = ","))
      ) +
    gb_theme
  
  ggsave(file.path(output_dir, paste0("Fig2_Extremes_", ds_name, ".png")), p2, width = 180, height = 120, units = "mm", dpi = 600)

}

# ------------------------------------------------------------------------------
# 5. MULTI-PANEL A4 GRID PLOTTING LOOP
# ------------------------------------------------------------------------------
plot_list <- list()

# Shared Theme for Sub-plots
sub_theme <- theme_bw(base_size = 6) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    panel.grid.major = element_line(color = "grey95"), 
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 6, color = "black"),
    
    # Tiny margin to fit in grid
    plot.margin = margin(2, 2, 2, 12, "mm"),
    
    # X Axis: Small rotated text
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color="black"),
    axis.text.y = element_text(size = 6, color="black"),
    axis.title = element_blank(), # Remove axis titles for cleaner grid
    
    legend.position = "none", # Hide legend in individual plots
    plot.title = element_text(face = "bold", size = 8, hjust = 0, margin = margin(b=2))
  )

for (ds_name in unique_datasets) {
  
  # Filter: Common + All Rare
  p_data <- df_final %>%
    filter(dataset == ds_name) %>%
    # filter(!is.na(abundance_class)) %>%
    mutate(
      abundance_class = factor(abundance_class, levels = c("Most prevalent", "Least prevalent")),
      pop_label = reorder(pop_label, -actual_count)
    )
  
  ds_name_n_cells <- df_final %>% filter(dataset == ds_name) %>% select(dataset_total) %>% unique() %>% pull()
  
  # Skip empty datasets
  if(nrow(p_data) == 0) next
  
  p <- ggplot(p_data, aes(x = pop_label, y = f1_score, fill = model, color = model)) +
    geom_boxplot(
      outlier.size = 0.05, 
      lwd = 0.25, 
      alpha = 0.7,
      width = 0.7,
      fatten = 0.8
    ) +
    # facet_grid(. ~ abundance_class, scales = "free_x", space = "free_x") +
    
    scale_fill_manual(values = tool_colors) +
    scale_color_manual(values = tool_colors) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    
    labs(title = ds_name, subtitle = paste("N cells:", format(ds_name_n_cells, big.mark = ","))) +
    sub_theme
  
  plot_list[[ds_name]] <- p
}

# ASSEMBLE GRID (3 LEFT, 3 RIGHT = 2 Columns)
plot_list <- plot_list[sort(names(plot_list))]

legend_plot <- ggplot(df_final, aes(x=model, y=f1_score, fill=model, color=model)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values = tool_colors, name = "Method") +
  scale_color_manual(values = tool_colors, name = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size=8, face="bold"))

common_legend <- cowplot::get_legend(legend_plot)

final_grid <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(
    title = "",
    theme = theme(plot.title = element_text(face="bold", size=12))
  )

final_figure <- cowplot::plot_grid(
  final_grid,
  common_legend,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

ggsave(
  file.path(output_dir, "Figure2_MultiPanel_A4.png"),
  plot = final_figure,
  width = 220,
  height = 300,
  units = "mm",
  dpi = 2000
)

print(paste("Saved to", file.path(output_dir, "Figure2_MultiPanel_A4.png")))



# # 1to4
# final_grid <- wrap_plots(plot_list[1:4], ncol = 1) +
#   plot_annotation(
#     title = "",
#     theme = theme(plot.title = element_text(face="bold", size=12))
#   )
# 
# final_figure <- cowplot::plot_grid(
#   final_grid,
#   common_legend,
#   ncol = 1,
#   rel_heights = c(1, 0.05)
# )
# 
# ggsave(
#   file.path(output_dir, "Figure2_MultiPanel_A4_1to4.png"),
#   plot = final_figure,
#   width = 220,
#   height = 300,
#   units = "mm",
#   dpi = 2000
# )
# 
# # 5to8
# final_grid <- wrap_plots(plot_list[5:8], ncol = 1) +
#   plot_annotation(
#     title = "",
#     theme = theme(plot.title = element_text(face="bold", size=12))
#   )
# 
# final_figure <- cowplot::plot_grid(
#   final_grid,
#   common_legend,
#   ncol = 1,
#   rel_heights = c(1, 0.05)
# )
# 
# ggsave(
#   file.path(output_dir, "Figure2_MultiPanel_A4_5to8.png"),
#   plot = final_figure,
#   width = 220,
#   height = 300,
#   units = "mm",
#   dpi = 2000
# )


# ------------------------------------------------------------------------------
# 6. REFINED A4 GRID (Top 5 Common & Rare)
# ------------------------------------------------------------------------------
# FILTER DATASETS (Must have both Common and Rare)
# valid_datasets <- pop_meta %>%
#   group_by(dataset) %>%
#   # summarize(
#   #   has_common = any(abundance_class == "Most prevalent"),
#   #   has_rare   = any(abundance_class == "Least prevalent")
#   # ) %>%
#   # filter(has_common & has_rare) %>%
#   pull(dataset)

# valid_datasets <- pop_meta$dataset %>% unique()
valid_datasets <- c("COVID", "CV", "HBM", "LV", "MBM", "SB")

# SUBPLOT THEME
gb_theme_subplot <- theme_bw(base_size = 8) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.2), 
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 7, color = "black"),
    
    # CRITICAL: Margins to prevent name cutting
    plot.margin = margin(2, 2, 2, 16, "mm"),
    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5, color="black"),
    axis.text.y = element_text(size = 6, color="black"),
    axis.title = element_blank(),
    
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 9, hjust = 0, margin = margin(b=3))
  )

plot_list_refined <- list()

for (ds_name in valid_datasets) {
  
  # ds_name <- "COVID"
  
  p_data <- df_final %>%
    filter(dataset == ds_name) %>%
    filter(!is.na(abundance_class)) %>%
    filter(rank_within_class <= 5) %>% 
    mutate(
      abundance_class = factor(abundance_class, levels = c("Most prevalent", "Least prevalent")),
      pop_label_newline = reorder(pop_label_newline, -pct)
    )
  
  ds_name_n_cells <- df_final %>% filter(dataset == ds_name) %>% select(dataset_total) %>% unique() %>% pull()
  
  n_pop <- p_data$population_name %>% unique() %>% length()
 
  p <- ggplot(p_data, aes(x = pop_label_newline, y = f1_score, fill = model, color = model)) +
    geom_boxplot(
      outlier.size = 0.1, 
      lwd = 0.25,      
      alpha = 0.7,     
      width = 0.65,    
      fatten = 1       
    ) +
    
    (if (n_pop >= 6) facet_grid(. ~ abundance_class, scales = "free_x", space = "free_x") else NULL) +
    
    # facet_grid(. ~ abundance_class, scales = "free_x", space = "free_x") +
    
    scale_fill_manual(values = tool_colors) +
    scale_color_manual(values = tool_colors) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    
    labs(title = ds_name, subtitle = paste("N cells:", format(ds_name_n_cells, big.mark = ","))) +
    gb_theme_subplot
  
  plot_list_refined[[ds_name]] <- p
}

if(length(plot_list_refined) > 0) {
  
  legend_plot <- ggplot(df_final, aes(x=model, y=f1_score, fill=model, color=model)) +
    geom_boxplot(alpha=0.7) +
    scale_fill_manual(values = tool_colors, name = "Method") +
    scale_color_manual(values = tool_colors, name = "Method") +
    theme_minimal() +
    theme(
      legend.position = "bottom", 
      legend.text = element_text(size=8, face="bold"),
      legend.key.size = unit(4, "mm")
    )
  
  common_legend <- get_legend(legend_plot)
  
  p_grid <- plot_grid(plotlist = plot_list_refined, ncol = 2, align = 'hv', axis = 'tb')
  
  final_plot <- plot_grid(
    p_grid,
    common_legend,
    ncol = 1,
    rel_heights = c(1.5, 0.05) 
  )
  
  outfile <- file.path(output_dir, "Figure2_A4_Top5_Refined.png")
  
  ggsave(
    outfile, 
    plot = final_plot, 
    # width = 210, 
    # height = 297, 
    width = 210, 
    height = 297, 
    units = "mm", 
    dpi = 2000
  )
  
   print(paste("Saved A4 Figure to:", outfile))
   

   # p_grid <- plot_grid(plotlist = plot_list_refined[1:4], ncol = 1, align = 'hv', axis = 'tb')
   # 
   # final_plot <- plot_grid(
   #   p_grid,
   #   common_legend,
   #   ncol = 1,
   #   rel_heights = c(1, 0.05) 
   # )
   # 
   # ggsave(
   #   file.path(output_dir, "Figure2_A4_Top5_Refined_1to4.png"), 
   #   plot = final_plot, 
   #   # width = 210, 
   #   # height = 297, 
   #   width = 210, 
   #   height = 297, 
   #   units = "mm", 
   #   dpi = 2000
   # )
   # 
   # plot_list_p2 <- plot_list_refined[5:length(plot_list_refined)]
  
} else {
  print("No valid datasets found.")
}

###
# Usage example:
#Rscript Fig2_plot.r \
#  --conf_input ../ob-blob-metrics/out/metric_collectors/metrics_report/per_population_confusion.tsv \
#  --output_dir ../ob-pipeline-plots
###