#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Main configuration file path
params.config_file = "./example_config.json"
params.outdir = "./results"
params.epiexpress_version = null

// Import pipeline modules
include { SETUP_EPIEXPRESS } from "./modules/setup_epiexpress.nf"
include { EPIEXPRESS_PREDICTION } from "./modules/epiexpress_prediction.nf"
include { GENERATE_MATRIX as GENERATE_MATRIX_RF } from "./modules/generate_matrix.nf"
include { GENERATE_MATRIX as GENERATE_MATRIX_CNN } from "./modules/generate_matrix.nf"
include { DIFFERENTIAL_EXPRESSION } from "./modules/differential_expression.nf"
include { REAL_RNA_COMPARISON } from "./modules/real_rna_comparison.nf"
include { GO_ENRICHMENT } from "./modules/go_enrichment.nf"
include { DISGENET } from "./modules/disgenet.nf"

workflow {
    // Load configuration from JSON file
    def config_file = file(params.config_file)
    def common_config = new groovy.json.JsonSlurper().parseText(config_file.text)

    // Use out_folder from config if specified, otherwise use command line param
    def base_output_dir = common_config.EpiExpress.out_folder
        ? file(common_config.EpiExpress.out_folder).toAbsolutePath().toString()
        : file(params.outdir).toAbsolutePath().toString()

    // Create output directory based on design formula
    def design_name = common_config.HAWAS.design_formula
        .replaceAll(/^~\s*/, '')
        .replaceAll(/\s*\+\s*/, '_')
        .replaceAll(/[^\w_]/, '')

    def output_dir = "${base_output_dir}/${design_name}"
    log.info("Output: ${output_dir}")

    // Setup EpiExpress repository
    def config_repo = common_config.EpiExpress.repo_dir

    if (config_repo && file(config_repo).exists()) {
        epiexpress_channel = Channel.of(config_repo)
    }
    else {
        epiexpress_channel = SETUP_EPIEXPRESS()
    }

    // Stage 0: EpiExpress gene expression prediction
    EPIEXPRESS_PREDICTION(
        epiexpress_channel,
        common_config.EpiExpress.bigwigs,
        common_config.EpiExpress.cre_rf_model_folder,
        common_config.EpiExpress.binned_cnn_model_folder,
        common_config.EpiExpress.provided_input,
        common_config.EpiExpress.gene_file,
        common_config.EpiExpress.mode,
        common_config.EpiExpress.correlation_cutoff,
        common_config.EpiExpress.cores,
    )

    // Stage 1: Generate matrices
    def mode = common_config.EpiExpress.mode
    def matrix_cores = common_config.EpiExpress.cores ?: 1
    rf_matrix = Channel.empty()
    cnn_matrix = Channel.empty()

    if (mode == 'CRE-RF' || mode == 'all' || mode == 'CRE-RF,Binned-CNN') {
        GENERATE_MATRIX_RF(
            EPIEXPRESS_PREDICTION.out.map { it.resolve("CRE_RF_result") },
            "${base_output_dir}/matrices",
            matrix_cores,
            "combined_matrix_CRE-RF.csv.gz",
        )
        rf_matrix = GENERATE_MATRIX_RF.out.matrix_file
    }

    if (mode == 'Binned-CNN' || mode == 'all' || mode == 'CRE-RF,Binned-CNN') {
        GENERATE_MATRIX_CNN(
            EPIEXPRESS_PREDICTION.out.map { it.resolve("Binned_CNN_Result") },
            "${base_output_dir}/matrices",
            matrix_cores,
            "combined_matrix_Binned-CNN.csv.gz",
        )
        cnn_matrix = GENERATE_MATRIX_CNN.out.matrix_file
    }

    // Stage 2: HAWAS Differential Expression Analysis
    DIFFERENTIAL_EXPRESSION(
        rf_matrix.ifEmpty([]),
        cnn_matrix.ifEmpty([]),
        output_dir,
        file(common_config.HAWAS.metadata_file),
        common_config.HAWAS.design_formula,
        common_config.HAWAS.metadata_sample_id,
    )

    // Stage 3 (optional): Real RNA comparison
    def real_rna_path = common_config.HAWAS.get('real_rna_matrix')
    if (real_rna_path == 'null') {
        real_rna_path = null
    }
    real_sig_genes = Channel.empty()
    if (real_rna_path) {
        REAL_RNA_COMPARISON(
            Channel.of(file(real_rna_path)),
            DIFFERENTIAL_EXPRESSION.out.sig_genes.filter { it.parent.name == "CRE-RF" }.ifEmpty([]),
            DIFFERENTIAL_EXPRESSION.out.sig_genes.filter { it.parent.name == "Binned-CNN" }.ifEmpty([]),
            rf_matrix.ifEmpty([]),
            cnn_matrix.ifEmpty([]),
            output_dir,
            file(common_config.HAWAS.metadata_file),
            common_config.HAWAS.design_formula,
            common_config.HAWAS.metadata_sample_id,
        )
        real_sig_genes = REAL_RNA_COMPARISON.out.sig_genes
    }

    // Stage 4 (optional): GO enrichment per model
    // Runs only if "run_go_enrichment": true is set in HAWAS config
    def gene_list_file = file(common_config.EpiExpress.gene_file)
    if (common_config.HAWAS.get('run_go_enrichment') == true) {
        def sig_genes_with_model = DIFFERENTIAL_EXPRESSION.out.sig_genes.map { sig_genes_file ->
            tuple(sig_genes_file, sig_genes_file.parent.name)
        }
        GO_ENRICHMENT(sig_genes_with_model, output_dir, gene_list_file)
    }


    // Stage 5 (optional): DISGENET disease enrichment
    // runs only if a GMT file is provided and all three gene lists are available
    // (CRE-RF, Binned-CNN, and Real-RNA).
    def disgenet_gmt = common_config.HAWAS.get('disgenet_gmt')
    def has_all_models = (mode == 'all' || mode == 'CRE-RF,Binned-CNN') && real_rna_path
    if (disgenet_gmt && has_all_models) {
        DISGENET(
            DIFFERENTIAL_EXPRESSION.out.sig_genes.filter { it.parent.name == "CRE-RF" },
            DIFFERENTIAL_EXPRESSION.out.sig_genes.filter { it.parent.name == "Binned-CNN" },
            real_sig_genes,
            gene_list_file,
            file(disgenet_gmt),
            output_dir,
        )
    }
    else if (disgenet_gmt) {
        log.warn("DISGENET skipped: requires mode='all' (or 'CRE-RF,Binned-CNN') and a real_rna_matrix.")
    }
}
