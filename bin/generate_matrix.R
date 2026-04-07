#!/usr/bin/env Rscript

library(argparse)
library(parallel)
library(data.table)

parser <- ArgumentParser()
parser$add_argument("--input_dir",
  required = TRUE,
  help = "Directory with per-gene prediction CSVs"
)
parser$add_argument("--output_file",
  required = TRUE,
  help = "Output filename (.csv.gz)"
)
parser$add_argument("--n_cores",
  type = "integer", default = 1,
  help = "Cores for mclapply"
)
args <- parser$parse_args()

target_column <- "Backscaled_prediction"

csv_files <- list.files(args$input_dir,
  pattern = "_predictions\\.csv\\.gz$", full.names = TRUE
)

if (length(csv_files) == 0) stop("No prediction files found in ", args$input_dir)
cat("Found", length(csv_files), "prediction files\n")

gene_names <- sub("_predictions\\.csv\\.gz$", "", basename(csv_files))

# sample names from first file
ref_dt <- fread(csv_files[1])
sample_ids <- ref_dt$Sample
cat("Samples:", length(sample_ids), "\n")
cat("Genes:", length(gene_names), "\n")

results <- mclapply(csv_files, function(f) {
  fread(f, select = target_column)[[target_column]]
}, mc.cores = args$n_cores)

final_matrix <- do.call(rbind, results)
rownames(final_matrix) <- gene_names
colnames(final_matrix) <- sample_ids

keep <- complete.cases(final_matrix)
cat("Removed", sum(!keep), "genes with NAs\n")
final_matrix <- final_matrix[keep, ]

output_dt <- as.data.table(final_matrix, keep.rownames = TRUE)
fwrite(output_dt, file = args$output_file, compress = "gzip")

cat("Matrix:", nrow(final_matrix), "genes x", ncol(final_matrix), "samples\n")
cat("Output:", args$output_file, "\n")
