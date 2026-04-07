nextflow.enable.dsl = 2

process REAL_RNA_COMPARISON {
    tag "Real_RNA_Comparison"
    publishDir "${output_dir}", mode: 'copy'

    input:
    path real_rna_matrix
    path rf_deseq2_results
    path cnn_deseq2_results
    path rf_matrix
    path cnn_matrix
    val output_dir
    path metadata_file
    val design_formula
    val metadata_sample_id

    output:
    path "DESeq2_results/Real-RNA"
    path "DESeq2_results/Comparison", optional: true
    path "DESeq2_results/Real-RNA/sig_genes.txt", emit: sig_genes

    script:
    def rf_res = rf_deseq2_results ? "--rf_deseq2_results ${rf_deseq2_results}" : ""
    def cnn_res = cnn_deseq2_results ? "--cnn_deseq2_results ${cnn_deseq2_results}" : ""
    def rf_mat = rf_matrix ? "--rf_matrix ${rf_matrix}" : ""
    def cnn_mat = cnn_matrix ? "--cnn_matrix ${cnn_matrix}" : ""
    """
    Rscript ${projectDir}/bin/real_rna_comparison.R \
        --real_rna_matrix ${real_rna_matrix} \
        --metadata_file ${metadata_file} \
        --design_formula "${design_formula}" \
        --metadata_sample_id ${metadata_sample_id} \
        --output_dir . \
        ${rf_res} ${cnn_res} ${rf_mat} ${cnn_mat}
    """
}
