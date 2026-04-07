#!/usr/bin/env Rscript

library(gprofiler2)
library(data.table)
library(ggplot2)

library(argparse)

parser <- ArgumentParser()
parser$add_argument("--background",
  required = TRUE,
  help = "Gene background list"
)
parser$add_argument("--deseq_results",
  required = TRUE,
  help = "DESeq2 results table"
)
args <- parser$parse_args()

background_file <- args$background
deseq_file <- args$deseq_results

background <- unique(fread(background_file, header = FALSE)[[1]])
deseq_res <- fread(deseq_file)

sig <- deseq_res[
  !is.na(padj) & !is.na(log2FoldChange) &
    padj <= 0.05 & abs(log2FoldChange) >= 0.585 & baseMean >= 10
]

if (nrow(sig) == 0) {
  message("No significant genes")
  quit(save = "no", status = 0)
}

message(sprintf("%d significant genes total", nrow(sig)))
message(sprintf("%d upregulated", sum(sig$log2FoldChange > 0)))
message(sprintf("%d downregulated", sum(sig$log2FoldChange < 0)))

setorder(sig, padj)
up_genes <- sig[log2FoldChange > 0, gene]
down_genes <- sig[log2FoldChange < 0, gene]

run_gost <- function(genes, name) {
  if (length(genes) < 3) {
    return()
  }

  res <- gost(genes,
    organism = "hsapiens", significant = TRUE,
    user_threshold = 0.05, correction_method = "g_SCS",
    custom_bg = background, domain_scope = "custom",
    sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
    ordered_query = TRUE, exclude_iea = TRUE
  )

  if (is.null(res$result) || nrow(res$result) == 0) {
    return()
  }

  fwrite(res$result, paste0(name, "_enrichment.txt"), sep = "\t", quote = FALSE)

  top <- as.data.table(res$result)
  setorder(top, p_value)
  top <- top[term_size >= 10 & term_size <= 500]
  top <- top[, .SD[seq_len(min(5, .N))], by = source]
  top[, term_name := strtrim(term_name, 45)]

  if (nrow(top) == 0) {
    return()
  }

  src_labels <- c(
    "GO:BP" = "GO: Biological Process",
    "GO:MF" = "GO: Molecular Function",
    "GO:CC" = "GO: Cellular Component",
    "KEGG" = "KEGG", "REAC" = "Reactome"
  )

  p <- ggplot(top, aes(x = -log10(p_value), y = term_name, fill = source)) +
    geom_col() +
    facet_wrap(~source, scales = "free_y", ncol = 1, labeller = labeller(source = src_labels)) +
    labs(x = "-log10(p-value)", y = "", title = paste(name, "regulated")) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 9), legend.position = "none")

  ggsave(paste0(name, "_enrichment.pdf"), p,
    width = 12,
    height = max(6, nrow(top) * 0.45)
  )
}

run_gost(up_genes, "up")
run_gost(down_genes, "down")
message("Done")
