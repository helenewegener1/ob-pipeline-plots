#!/usr/bin/env bash
set -e

echo "Generating Figure 1..."
## Generate Figure 1

Rscript Fig1_plot.r \
  --f1_input ../ob-blob-metrics/out/metric_collectors/metrics_report/f1_macro_by_crossvalidation.tsv \
  --output ../ob-pipeline-plots/Figure1_Heatmap-boxplot.png

echo "Generating Figure 2..."
## Generate Figure 2

Rscript Fig2_plot.r \
  --conf_input ../ob-blob-metrics/out/metric_collectors/metrics_report/per_population_confusion.tsv \
  --output_dir ../ob-pipeline-plots

echo "Generating Figure 3..."
## Generate Figure 3

Rscript combined_plot_3.r \
   --f1_weighted ../ob-blob-metrics/out/metric_collectors/metrics_report/samples_vs_f1_weighted.tsv \
   --perf_input ../ob-blob-metrics/out/performances.tsv \
   --meta_json ../ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
   --output ../ob-pipeline-plots/Figure3_Combined_Grid.png

echo "All plots generated successfully."

## Example usage:
