nextflow.enable.dsl = 2

process DISGENET {
    tag "DISGENET"
    publishDir "${output_dir}/disgenet", mode: 'copy'

    input:
    path rf_sig_genes
    path cnn_sig_genes
    path real_sig_genes
    path background
    path gmt_file
    val output_dir

    output:
    path "*.pdf", optional: true

    script:
    """
    Rscript ${projectDir}/bin/disgenet.R \
        --gmt_file ${gmt_file} \
        --background ${background} \
        --rf_genes ${rf_sig_genes} \
        --cnn_genes ${cnn_sig_genes} \
        --real_genes ${real_sig_genes}
    """
}
