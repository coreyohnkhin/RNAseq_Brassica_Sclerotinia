#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(msigdbr)
  library(fgsea); library(biomaRt)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Pass one or more results_*.csv files")

mart_bo <- useMart("plants_mart", dataset="bo_oleracea_eg_gene")
mart_at <- useMart("plants_mart", dataset="athaliana_eg_gene")

map <- getLDS(attributes="ensembl_gene_id", mart=mart_bo,
              attributesL="ensembl_gene_id", martL=mart_at,
              uniqueRows=TRUE) |>
       rename(bo_gene = Gene.stable.ID, at_gene = Gene.stable.ID.1)

ats_go <- msigdbr(species="Arabidopsis thaliana", category="C5") |>
          split(.$entrez_gene)

make_ranks <- function(csv, stat_col="stat") {
  fread(csv) |> arrange(desc(!!sym(stat_col))) |>
    filter(!is.na(!!sym(stat_col))) |> {
      setNames(.$[[stat_col]], .$gene)
    }
}

for (f in args) {
  message("▶ Processing ", f)
  ranks <- make_ranks(f, "stat")

  ranks_at <- tibble(gene=names(ranks), stat=unname(ranks)) |>
              inner_join(map, by=c(gene="bo_gene")) |>
              transmute(at_gene, stat) |> deframe()

  fg <- fgsea(pathways=ats_go, stats=ranks_at,
              nperm=10000, minSize=15, maxSize=500, eps=0.0) |>
        arrange(padj)

  out <- file.path(dirname(f), paste0("fgsea_", basename(f) |> tools::file_path_sans_ext(), ".tsv"))
  fwrite(fg, out, sep="\t")
  message("    → written ", out)
}
message("All done ✓")
