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
  library(patchwork)
  library(stringr)
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-f", "--f1_input"), type = "character", default = NULL, 
              help = "Path to samples_vs_f1_weighted.tsv", metavar = "character"),
  make_option(c("-p", "--perf_input"), type = "character", default = NULL, 
              help = "Path to performances.tsv", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "Figure5_Scalability.png", 
              help = "Output filename", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validation
if (is.null(opt$f1_input) || is.null(opt$perf_input)) {
  print_help(opt_parser)
  stop("Both --f1_input and --perf_input must be provided.", call. = FALSE)
}

path_f1 <- opt$f1_input
path_perf <- opt$perf_input
output_dir <- opt$output

# ------------------------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------------------------
df_f1_raw <- read_tsv(path_f1, show_col_types = FALSE)
df_perf_raw <- read_tsv(path_perf, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 2. DATA PROCESSING (LEVINE ONLY)
# ------------------------------------------------------------------------------

# --- MAPPINGS (Consistent with Figs 2, 3, 4) ---
model_map <- c(
  'cyanno'      = "CyAnno",
  'cygate'      = "CyGATE",
  'dgcytof'     = "DGCytof",
  'gatemeclass' = "GateMeClass",
  'lda'         = "CyTOF LC",
  'random'      = "Random"
)

# --- COLORS (Set1 Palette) ---
tool_colors <- c(
  "CyAnno"      = "#E41A1C",  # Red
  "CyGATE"      = "#377EB8",  # Blue
  "DGCytof"     = "#4DAF4A",  # Green
  "GateMeClass" = "#984EA3",  # Purple
  "CyTOF LC"    = "#FF7F00",  # Orange
  "Random"      = "#525252"   # Dark Grey
)

# 2. DATA PROCESSING (Using your robust logic + Adding Name Mapping)

# Helper to extract size
extract_size <- function(name) {
  as.numeric(str_extract(name, "(?<=sub-sampling-)\\d+"))
}

# A. Process Performance (Time)
df_levine_time <- df_perf_raw %>%
  filter(str_detect(params, "sub-sampling")) %>%
  mutate(
    # Your extraction logic
    train_size = as.numeric(str_match(params, 'sub-sampling[^:]+:\\s*\"+([^\"]+)')[,2]),
    # Clean Tool Name
    model = recode(module, !!!model_map)
  ) %>%
  filter(!module %in% c("metrics", "flow_metrics", "data_import", "data_preprocessing", "random")) %>%
  filter(model %in% names(tool_colors)) %>% # Safety filter
  group_by(model, train_size) %>%
  summarise(mean_time = mean(s, na.rm = TRUE), .groups = "drop")

# B. Process F1 (Accuracy)
df_levine_f1 <- df_f1_raw %>%
  filter(str_detect(dataset, "sub-sampling")) %>%
  mutate(
    # Your extraction logic
    train_size = extract_size(dataset),
    # Clean Tool Name
    model = recode(model, !!!model_map)
  ) %>%
  filter(!str_detect(model, "Random")) %>%
  filter(model %in% names(tool_colors)) %>%
  group_by(model, train_size) %>%
  summarise(mean_f1 = mean(f1_weighted_median, na.rm = TRUE), .groups = "drop")


# 3. PLOTTING SETUP (Visual Identity from Figs 3 & 4)

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
p_time <- ggplot(df_levine_time, aes(x = train_size, y = mean_time, color = model, group = model)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  # Points: Filled with color, thin black border (Matches Fig 3 Scatter style)
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
p_f1 <- ggplot(df_levine_f1, aes(x = train_size, y = mean_f1, color = model, group = model)) +
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
final_fig5 <- (p_time / p_f1) + 
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

outfile <- file.path(output_dir, "Figure5_Healthy_Scalability_Final.png")

ggsave(outfile, 
       plot = final_fig5, 
       width = 120, 
       height = 150, 
       units = "mm", 
       dpi = 800)

print(paste("Saved Figure 5 to", outfile))

###
# Usage example:
#Rscript Fig5_plot.r \
#  --f1_input ./ob-blob-metrics/out/metric_collectors/metrics_report/samples_vs_f1_weighted.tsv \
#  --perf_input ./ob-blob-metrics/out/performances.tsv \
#  --output ./ob-pipeline-plots/
###