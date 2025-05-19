#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("DESeq2","ashr","tidyverse","ComplexHeatmap")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

library(DESeq2); library(ashr); library(tidyverse); library(ComplexHeatmap)

ctrl <- read_delim("data/featurecountmatrices/N1909Pet_ControlBOL.txt",
                   skip = 1, show_col_types = FALSE)
inf  <- read_delim("data/featurecountmatrices/N1909Pet_InfectedBOL.txt",
                   skip = 1, show_col_types = FALSE)

drop_annot <- c("Chr","Start","End","Strand","Length")
counts <- inner_join(ctrl  |> select(-all_of(drop_annot)),
                     inf   |> select(-all_of(drop_annot)),
                     by = "Geneid")

dds <- DESeqDataSetFromMatrix(countData = as.matrix(column_to_rownames(counts,"Geneid")),
                              colData   = tibble(sample = colnames(counts[-1]),
                                                 condition = rep(c("control","infected"),
                                                                 c(ncol(ctrl)-1,ncol(inf)-1)))
                                           |> column_to_rownames("sample"),
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) > 10, ]
res <- lfcShrink(DESeq(dds), coef="condition_infected_vs_control", type="ashr") |>
       as_tibble(rownames="gene") |> arrange(padj)

dir.create("results", showWarnings = FALSE)
write_csv(res,"results/results_BOL_infected_vs_control.csv")
cat("Significant genes (padj < 0.05):",sum(res$padj < 0.05,na.rm=TRUE),"\n")

vst <- varianceStabilizingTransformation(dds) |> assay()
Heatmap(vst[head(res$gene,50),], name="VST",
        column_split = dds$condition,
        top_annotation = HeatmapAnnotation(condition = dds$condition)) |>
  (function(ht){
     dir.create("figures", showWarnings = FALSE)
     pdf("figures/heatmap_top50_BOL.pdf",6,8); draw(ht); dev.off()
   })
