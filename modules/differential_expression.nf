nextflow.enable.dsl = 2

process DIFFERENTIAL_EXPRESSION {
    tag "HAWAS_Differential_Expression"
    publishDir "${output_dir_deseq}", mode: 'copy'

    input:
    path combined_matrix_RF
    path combined_matrix_CNN
    val output_dir_deseq
    path metadata_file
    val design_formula
    val metadata_sample_id

    output:
    path "DESeq2_results"
    path "DESeq2_results/*/sig_genes.txt", emit: sig_genes

    script:
    def rf_matrix = combined_matrix_RF ? "${combined_matrix_RF}" : ""
    def cnn_matrix = combined_matrix_CNN ? "${combined_matrix_CNN}" : ""
    """

    Rscript ${projectDir}/bin/deseq2_analysis.R \
        --combined_matrix_rf "${rf_matrix}" \
        --combined_matrix_cnn "${cnn_matrix}" \
        --output_dir . \
        --metadata_file "${metadata_file}" \
        --design_formula "${design_formula}" \
        --metadata_sample_id "${metadata_sample_id}"
    """
}
