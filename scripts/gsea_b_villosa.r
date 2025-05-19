#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  cat("\nGene‑set enrichment analysis (fgsea)\n",
      "USAGE: Rscript scripts/gsea.r <results_csv1> [<results_csv2> …]\n",
      "Each <results_csv> must be a DESeq2 results table with a ‘gene’ column\n",
      "and Wald statistic ‘stat’.\n\n", sep = "")
  quit(status = 1)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

pkgs <- c("fgsea", "data.table", "tidyverse", "msigdbr", "biomaRt")
for (p in pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(fgsea); library(data.table); library(tidyverse)
  library(msigdbr); library(biomaRt)
})

make_ranks <- function(csv, rank_col = c("stat", "log2FoldChange", "signed_logp")) {
  rank_col <- match.arg(rank_col)
  dt <- fread(csv) |> drop_na(pvalue)
  if (rank_col == "signed_logp")
    dt <- mutate(dt, signed_logp = sign(log2FoldChange) * -log10(pvalue))
  deframe(dt |> arrange(desc(.data[[rank_col]])) |> select(gene, .data[[rank_col]]))
}

message("Connecting to Ensembl Plants (biomaRt)…")
mart_bo <- useMart("plants_mart", dataset = "bo_oleracea_eg_gene")
mart_at <- useMart("plants_mart", dataset = "athaliana_eg_gene")

message("Retrieving orthology map…")
map <- getLDS(attributes = "ensembl_gene_id",        mart  = mart_bo,
              attributesL = "ensembl_gene_id",       martL = mart_at,
              uniqueRows  = TRUE) |>
       rename(bo_gene = Gene.stable.ID, at_gene = Gene.stable.ID.1)

message("Downloading Arabidopsis GO gene sets (MSigDB)…")
ats_go <- msigdbr(species = "Arabidopsis thaliana", category = "C5") |>
          split(.$entrez_gene)

run_fgsea <- function(res_csv) {
  stopifnot(file.exists(res_csv))
  message("▶  Processing ", res_csv)

  ranks <- make_ranks(res_csv, "stat")

  ranks_at <- tibble(gene = names(ranks), stat = unname(ranks)) |>
              inner_join(map, by = c(gene = "bo_gene")) |>
              transmute(at_gene, stat) |>
              deframe()

  fg <- fgsea(pathways = ats_go,
              stats    = ranks_at,
              nperm    = 10000,
              minSize  = 15,
              maxSize  = 500,
              eps      = 0.0) |>
        arrange(padj)

  in_dir   <- dirname(res_csv)
  out_dir  <- ifelse(dir.exists(in_dir), in_dir, "results")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  contrast <- sub("^results_", "", tools::file_path_sans_ext(basename(res_csv)))
  out_file <- file.path(out_dir, paste0("fgsea_", contrast, ".tsv"))
  fwrite(fg, out_file, sep = "\t")
  message("    → written ", out_file)
}

invisible(lapply(args, run_fgsea))
message("All done ✓")
