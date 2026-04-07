library(argparse)
library(DESeq2)
library(data.table)
library(ggplot2)

parser <- ArgumentParser(description = "Real RNA DESeq2 + comparison with model predictions")
parser$add_argument("--real_rna_matrix", required = TRUE)
parser$add_argument("--rf_deseq2_results", help = "sig_genes.txt from CRE-RF (optional)")
parser$add_argument("--cnn_deseq2_results", help = "sig_genes.txt from Binned-CNN (optional)")
parser$add_argument("--output_dir", required = TRUE)
parser$add_argument("--metadata_file", required = TRUE)
parser$add_argument("--design_formula", required = TRUE)
parser$add_argument("--metadata_sample_id", required = TRUE)
parser$add_argument("--rf_matrix", help = "CRE-RF combined matrix CSV.GZ for scatter plots")
parser$add_argument("--cnn_matrix", help = "Binned-CNN combined matrix CSV.GZ for scatter plots")
args <- parser$parse_args()

# metadata
metadata <- read.table(args$metadata_file,
  header = TRUE, sep = "\t",
  quote = "", stringsAsFactors = FALSE
)
id_col <- args$metadata_sample_id
if (!id_col %in% colnames(metadata)) stop(sprintf("Column '%s' not in metadata.", id_col))
rownames(metadata) <- metadata[[id_col]]

out_base <- file.path(args$output_dir, "DESeq2_results")
rna_dir <- file.path(out_base, "Real-RNA")
if (!dir.exists(rna_dir)) dir.create(rna_dir, recursive = TRUE)

# load and align real RNA matrix
dt <- data.table::fread(args$real_rna_matrix, header = TRUE)
count_matrix <- as.matrix(dt[, -1, with = FALSE])
rownames(count_matrix) <- dt[[1]]

common_samples <- intersect(colnames(count_matrix), rownames(metadata))
if (length(common_samples) == 0) stop("No overlapping samples between RNA matrix and metadata.")
message(sprintf("Real RNA: %d samples matched to metadata.", length(common_samples)))

count_matrix <- count_matrix[, common_samples]
meta <- metadata[common_samples, ]
count_matrix <- count_matrix[, match(rownames(meta), colnames(count_matrix))]
stopifnot(all(rownames(meta) == colnames(count_matrix)))

# remove duplicate genes before DESeq2 renames them
count_matrix <- count_matrix[!duplicated(rownames(count_matrix)), ]

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData   = meta,
  design    = as.formula(args$design_formula)
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
message(sprintf("Real RNA: %d genes pass pre-filtering (rowSums >= 10).", sum(keep)))

# normalization for external raw RNA-seq counts (e.g. heart failure data)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

# full results + sig genes subset
res_df <- as.data.frame(res[order(res$padj), ])
res_df <- cbind(gene = rownames(res_df), res_df)
write.table(res_df, file.path(rna_dir, "deseq2_results.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

sig_df <- res_df[!is.na(res_df$padj) & res_df$padj <= 0.05 &
  abs(res_df$log2FoldChange) >= 0.585 & res_df$baseMean >= 10, ]
write.table(sig_df, file.path(rna_dir, "sig_genes.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message(sprintf("[Real-RNA] %d significant genes", nrow(sig_df)))

# plots
pdf(file.path(rna_dir, "dispersion_plot.pdf"), width = 8, height = 6)
plotDispEsts(dds)
dev.off()

pdf(file.path(rna_dir, "MA_plot.pdf"), width = 8, height = 6)
plotMA(res, ylim = c(-2, 2))
dev.off()

vsd <- vst(dds, blind = FALSE)
design_vars <- all.vars(as.formula(args$design_formula))
pdf(file.path(rna_dir, "PCA_plot.pdf"), width = 8, height = 6)
print(plotPCA(vsd, intgroup = design_vars[length(design_vars)]))
dev.off()

# scatter plots real RNA log2FC vs predicted log2FC
# sig genes from real RNA are the filter, log2FC computed from group means
# (same formula for both axes, Shamim's approach)
real_df <- as.data.frame(res)
sig_mask <- !is.na(real_df$padj) & real_df$padj <= 0.05 & abs(real_df$log2FoldChange) >= 0.585 & real_df$baseMean >= 10
sig_genes <- rownames(real_df)[sig_mask]

condition_col <- all.vars(as.formula(args$design_formula))
condition_col <- condition_col[length(condition_col)]
lvls <- levels(factor(meta[[condition_col]]))
ref_samp <- rownames(meta)[meta[[condition_col]] == lvls[1]]
treat_samp <- rownames(meta)[meta[[condition_col]] == lvls[2]]

rna_ref <- count_matrix[sig_genes, intersect(ref_samp, colnames(count_matrix)), drop = FALSE]
rna_treat <- count_matrix[sig_genes, intersect(treat_samp, colnames(count_matrix)), drop = FALSE]
log2fc_real_all <- log2((rowMeans(rna_ref) + 1) / (rowMeans(rna_treat) + 1))

model_matrices <- list()
if (!is.null(args$rf_matrix) && nzchar(args$rf_matrix) && file.exists(args$rf_matrix)) {
  model_matrices[["CRE-RF"]] <- args$rf_matrix
}
if (!is.null(args$cnn_matrix) && nzchar(args$cnn_matrix) && file.exists(args$cnn_matrix)) {
  model_matrices[["Binned-CNN"]] <- args$cnn_matrix
}

if (length(model_matrices) == 0 || length(sig_genes) == 0) {
  if (length(sig_genes) == 0) {
    message("No significant genes from real RNA, skipping scatter plots.")
  } else {
    message("No predicted matrices provided, skipping scatter plots.")
  }
} else {
  comp_dir <- file.path(out_base, "Comparison")
  if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)
  writeLines(sig_genes, file.path(comp_dir, "sig_genes_real.txt"))
  message(sprintf("Comparison: %d sig genes from real RNA.", length(sig_genes)))

  for (model_name in names(model_matrices)) {
    dt_model <- data.table::fread(model_matrices[[model_name]], header = TRUE)
    model_mat <- as.matrix(dt_model[, -1, with = FALSE])
    rownames(model_mat) <- dt_model[[1]]

    common_genes <- intersect(sig_genes, rownames(model_mat))
    if (length(common_genes) < 2) {
      message(sprintf("Skipping %s fewer than 2 sig genes in model matrix.", model_name))
      next
    }

    m_ref <- model_mat[common_genes, intersect(ref_samp, colnames(model_mat)), drop = FALSE]
    m_treat <- model_mat[common_genes, intersect(treat_samp, colnames(model_mat)), drop = FALSE]
    m_ref_means <- rowMeans(m_ref)
    m_treat_means <- rowMeans(m_treat)
    m_ref_means[m_ref_means == -1] <- 0 # -1 = failed prediction sentinel
    m_treat_means[m_treat_means == -1] <- 0
    log2fc_model <- log2((m_ref_means + 1) / (m_treat_means + 1))

    plot_data <- na.omit(data.frame(
      Gene         = common_genes,
      log2FC_real  = log2fc_real_all[common_genes],
      log2FC_model = log2fc_model
    ))
    if (nrow(plot_data) < 2) next

    r <- round(cor(plot_data$log2FC_real, plot_data$log2FC_model), 3)
    message(sprintf("Correlation Real vs %s: r=%.3f (n=%d genes)", model_name, r, nrow(plot_data)))

    p <- ggplot(plot_data, aes(x = log2FC_real, y = log2FC_model)) +
      geom_point(alpha = 0.4, size = 0.8) +
      xlim(-10, 10) +
      ylim(-10, 10) +
      labs(
        title = sprintf("Real RNA vs %s\n%d sig genes  |  r = %.3f", model_name, nrow(plot_data), r),
        x = "log2 Fold Change (Real RNA)",
        y = sprintf("log2 Fold Change (%s)", model_name)
      ) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
      )

    out_file <- sprintf("scatter_real_vs_%s.pdf", gsub("[^A-Za-z0-9]", "_", model_name))
    ggsave(file.path(comp_dir, out_file), plot = p, width = 6, height = 6)
  }
}

cat("\n=== Real RNA Comparison Complete ===\n")
