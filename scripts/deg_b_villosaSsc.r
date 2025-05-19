#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("DESeq2", "ashr", "tidyverse", "ComplexHeatmap")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}
library(DESeq2)
library(ashr)
library(tidyverse)
library(ComplexHeatmap)

ctrl_ssc <- read_delim(
  "data/featurecountmatrices/N1896Pet_ControlSsc.txt",
  skip   = 1,
  show_col_types = FALSE
)

inf_ssc  <- read_delim(
  "data/featurecountmatrices/N1896Pet_InfectedSsc.txt",
  skip   = 1,
  show_col_types = FALSE
)

drop_annot <- c("Chr", "Start", "End", "Strand", "Length")
counts_ctrl <- ctrl_ssc |> select(-all_of(drop_annot))
counts_inf  <- inf_ssc  |> select(-all_of(drop_annot))

counts <- inner_join(counts_ctrl, counts_inf, by = "Geneid")

cnt_mat <- counts |>
  column_to_rownames("Geneid") |>
  as.matrix()

sample_info <- tibble(
  sample    = colnames(cnt_mat),
  condition = rep(c("control", "infected"),
                  c(ncol(counts_ctrl) - 1, ncol(counts_inf) - 1))
) |>
  column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(
  countData = cnt_mat,
  colData   = sample_info,
  design    = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds   <- DESeq(dds)
res   <- lfcShrink(dds,
                   coef = "condition_infected_vs_control",
                   type = "ashr") |>
         as_tibble(rownames = "gene") |>
         arrange(padj)

dir.create("results",  showWarnings = FALSE)
write_csv(res, "results/results_Ssc_infected_vs_control.csv")

cat("Significant genes (padj < 0.05):",
    sum(res$padj < 0.05, na.rm = TRUE), "\n")

vst_mat <- varianceStabilizingTransformation(dds) |> assay()

top50 <- res |>
  filter(!is.na(padj)) |>
  slice_head(n = 50) |>
  pull(gene)

ha <- HeatmapAnnotation(condition = sample_info$condition)

ht <- Heatmap(vst_mat[top50, ],
              name                      = "VST",
              show_row_names            = FALSE,
              column_split              = sample_info$condition,
              top_annotation            = ha,
              clustering_distance_rows  = "euclidean",
              clustering_method_rows    = "complete")

dir.create("figures", showWarnings = FALSE)
pdf("figures/heatmap_top50_Ssc.pdf", width = 6, height = 8)
draw(ht)
dev.off()

# Outputs:
# results/results_Ssc_infected_vs_control.csv
# figures/heatmap_top50_Ssc.pdf
