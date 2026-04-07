#!/usr/bin/env Rscript
#-----------code with ordering disease and picking 5 category with higest enrichment (average among three methods)
# Load libraries
library(argparse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

parser <- ArgumentParser(description = "DISGENET disease enrichment")
parser$add_argument("--gmt_file", required = TRUE, help = "DISGENET GMT file")
parser$add_argument("--background", required = TRUE, help = "Background gene list")
parser$add_argument("--rf_genes", required = TRUE, help = "CRE-RF sig_genes.txt")
parser$add_argument("--cnn_genes", required = TRUE, help = "Binned-CNN sig_genes.txt")
parser$add_argument("--real_genes", required = TRUE, help = "Real-RNA sig_genes.txt")
args <- parser$parse_args()

# Read DISGENET GMT file
lines <- readLines(args$gmt_file)
result <- list()
for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    result[[length(result) + 1]] <- data.frame(Disease = parts[2], Gene = parts[-(1:2)])
}
AML_datasets <- do.call(rbind, result)

# Read gene lists
CNN_genes <- fread(args$cnn_genes)$gene
RF_genes <- fread(args$rf_genes)$gene
real_genes <- fread(args$real_genes)$gene
all_genes <- readLines(args$background)

# Calculate the unique sets for CNN and RF
common_CNN_RF <- intersect(RF_genes, CNN_genes)
unique_RF <- setdiff(RF_genes, union(real_genes, CNN_genes))
unique_CNN <- setdiff(CNN_genes, union(real_genes, RF_genes))

# Initialize enrichment_data to store results
enrichment_data <- data.frame(
    Disease = character(),
    # unique_RF_adj = numeric(),
    # unique_CNN_adj = numeric(),
    real_adj = numeric(),
    RF_adj = numeric(),
    CNN_adj = numeric(),
    stringsAsFactors = FALSE
)

# Define a function to calculate enrichment using Fisher's test
calc_enrichment <- function(gene_list, disease_genes, all_genes) {
    NotSigGenes <- setdiff(all_genes, gene_list)
    in_Disease_and_SigGenes <- length(intersect(disease_genes, gene_list))
    in_Disease_and_NotSigGenes <- length(intersect(disease_genes, NotSigGenes))
    in_NotDisease_and_SigGenes <- length(intersect(setdiff(all_genes, disease_genes), gene_list))
    in_NotDisease_and_NotSigGenes <- length(intersect(setdiff(all_genes, disease_genes), NotSigGenes))

    fisher_result <- fisher.test(
        matrix(
            c(
                in_Disease_and_SigGenes, in_Disease_and_NotSigGenes,
                in_NotDisease_and_SigGenes, in_NotDisease_and_NotSigGenes
            ),
            nrow = 2, byrow = TRUE
        )
    )

    return(p.adjust(fisher_result$p.value, method = "BH"))
}

top_diseases <- unique(AML_datasets$Disease)
AML_top_datasets <- AML_datasets

# Loop through each disease to calculate enrichment
for (disease in top_diseases) {
    disease_genes <- AML_top_datasets %>%
        filter(Disease == disease) %>%
        pull(Gene)

    # Calculate adjusted p-values for each gene list
    real_adj_pval <- calc_enrichment(real_genes, disease_genes, all_genes)
    # unique_RF_adj_pval <- calc_enrichment(unique_RF, disease_genes, all_genes)
    # unique_CNN_adj_pval <- calc_enrichment(unique_CNN, disease_genes, all_genes)
    RF_p_adj_value <- calc_enrichment(RF_genes, disease_genes, all_genes)
    CNN_p_adj_value <- calc_enrichment(CNN_genes, disease_genes, all_genes)

    # Store the -log10(adjusted p-values) in the enrichment_data dataframe
    enrichment_data <- rbind(enrichment_data, data.frame(
        Disease = disease,
        # unique_RF_adj = -log10(unique_RF_adj_pval + 1e-10),
        # unique_CNN_adj = -log10(unique_CNN_adj_pval + 1e-10),
        real_adj = -log10(real_adj_pval + 1e-10),
        RF_adj = -log10(RF_p_adj_value + 1e-10),
        CNN_adj = -log10(CNN_p_adj_value + 1e-10)
    ))
}

# Reshape the enrichment_data for plotting
enrichment_data_melt <- enrichment_data %>%
    pivot_longer(
        cols = c(
            # "unique_RF_adj", "unique_CNN_adj",
            "real_adj", "RF_adj", "CNN_adj"
        ),
        names_to = "List", values_to = "Enrichment"
    )

# Specify the desired order of the List factor
enrichment_data_melt$List <- factor(enrichment_data_melt$List,
    levels = c(
        # "unique_RF_adj", "unique_CNN_adj",
        "real_adj", "RF_adj", "CNN_adj"
    )
)

# Filter out rows with near-zero enrichment values
enrichment_data_melt <- enrichment_data_melt %>% filter(Enrichment > 1e-10)

A <- as.data.frame(enrichment_data_melt)
# #max
# A$Disease <- factor(A$Disease,
#                     levels = A %>%
#                       group_by(Disease) %>%
#                       mutate(MaxEnrichment = max(Enrichment, na.rm = TRUE)) %>%
#                       distinct(Disease, MaxEnrichment) %>%
#                       arrange(desc(MaxEnrichment)) %>%
#                       dplyr::pull(Disease) %>%
#                       rev())  # Reverse the order

# average
A$Disease <- factor(A$Disease,
    levels = A %>%
        group_by(Disease) %>%
        mutate(AvgEnrichment = mean(Enrichment, na.rm = TRUE)) %>% # Use mean instead of max
        distinct(Disease, AvgEnrichment) %>%
        arrange(desc(AvgEnrichment)) %>%
        dplyr::pull(Disease) %>%
        rev()
) # Reverse the order

A <- A[A$Disease %in% tail(levels(A$Disease), 20), ]

# Now plot again
ggplot(A, aes(x = List, y = Disease, size = Enrichment, color = Enrichment)) +
    geom_point() +
    scale_size_continuous(range = c(1, 10)) + # Adjust dot sizes
    scale_color_gradient(low = "lightgreen", high = "darkgreen") + # Gradient color
    theme_minimal() +
    labs(
        title = "Enrichment of Gene Lists (adj)",
        x = "Gene List",
        y = "Disease Dataset",
        size = "-log10(Adjusted adj_p-value)",
        color = "-log10(Adjusted adj_p-value)"
    ) +
    theme(
        axis.text.y = element_text(size = 8), # Adjust text size
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1)
    ) # Rotate x-axis labels

ggsave("disgenet_enrichment_full_top_20.pdf", width = 10, height = 8)
