library(argparse)
library(DESeq2)
library(data.table)

parser <- ArgumentParser(description = "DESeq2 on EpiExpress predicted expression matrices")
parser$add_argument("--combined_matrix_rf", help = "CRE-RF combined matrix CSV.GZ")
parser$add_argument("--combined_matrix_cnn", help = "Binned-CNN combined matrix CSV.GZ")
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--metadata_file", required = TRUE)
parser$add_argument("--design_formula", required = TRUE)
parser$add_argument("--metadata_sample_id", required = TRUE)
args <- parser$parse_args()

metadata <- read.table(args$metadata_file,
  header = TRUE, sep = "\t",
  quote = "", stringsAsFactors = FALSE
)
id_col <- args$metadata_sample_id
if (!id_col %in% colnames(metadata)) stop(sprintf("Column '%s' not in metadata.", id_col))
rownames(metadata) <- metadata[[id_col]]

out_base <- file.path(args$output_dir, "DESeq2_results")
if (!dir.exists(out_base)) dir.create(out_base, recursive = TRUE)

run_deseq <- function(matrix_path, model_name) {
  out_dir <- file.path(out_base, model_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  dt <- data.table::fread(matrix_path, header = TRUE)
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt[[1]]

  common <- intersect(colnames(mat), rownames(metadata))
  mat <- mat[, common]
  meta <- metadata[common, ]
  mat <- mat[, match(rownames(meta), colnames(mat))]
  stopifnot(all(rownames(meta) == colnames(mat)))

  # drop genes with values too large for integer conversion
  keep <- apply(mat, 1, max) <= .Machine$integer.max
  if (any(!keep)) message(sprintf("[%s] dropped %d genes exceeding integer max", model_name, sum(!keep)))
  mat <- mat[keep, ]

  dds <- DESeqDataSetFromMatrix(
    countData = round(mat),
    colData   = meta,
    design    = as.formula(args$design_formula)
  )

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  message(sprintf("[%s] %d genes after pre-filtering", model_name, sum(keep)))

  sizeFactors(dds) <- 1
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)

  # full results + sig genes subset 
  res_df <- as.data.frame(res[order(res$padj), ])
  res_df <- cbind(gene = rownames(res_df), res_df)
  write.table(res_df, file.path(out_dir, "deseq2_results.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  sig_df <- res_df[!is.na(res_df$padj) & res_df$padj <= 0.05 &
    abs(res_df$log2FoldChange) >= 0.585 & res_df$baseMean >= 10, ]
  write.table(sig_df, file.path(out_dir, "sig_genes.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  message(sprintf("[%s] %d significant genes", model_name, nrow(sig_df)))

  # plots
  pdf(file.path(out_dir, "dispersion_plot.pdf"), width = 8, height = 6)
  plotDispEsts(dds)
  dev.off()

  pdf(file.path(out_dir, "MA_plot.pdf"), width = 8, height = 6)
  plotMA(res, ylim = c(-2, 2))
  dev.off()

  design_vars <- all.vars(as.formula(args$design_formula))
  vsd <- vst(dds, blind = FALSE)
  pdf(file.path(out_dir, "PCA_plot.pdf"), width = 8, height = 6)
  print(plotPCA(vsd, intgroup = design_vars[length(design_vars)]))
  dev.off()
}

if (!is.null(args$combined_matrix_rf) && file.exists(args$combined_matrix_rf)) {
  message("Running DESeq2 on CRE-RF predictions...")
  run_deseq(args$combined_matrix_rf, "CRE-RF")
}

if (!is.null(args$combined_matrix_cnn) && file.exists(args$combined_matrix_cnn)) {
  message("Running DESeq2 on Binned-CNN predictions...")
  run_deseq(args$combined_matrix_cnn, "Binned-CNN")
}

cat("\n=== DESeq2 Analysis Complete ===\n")
