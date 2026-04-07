nextflow.enable.dsl = 2

process GO_ENRICHMENT {
  tag "${model}"
  publishDir "${output_dir}/GO_enrichment/${model}", mode: 'copy'

  input:
  tuple path(deseq_results), val(model)
  val output_dir
  path gene_list

  output:
  path "*.txt", optional: true
  path "*.pdf", optional: true

  script:
  """
  Rscript ${projectDir}/bin/go_enrichment_analysis.R \
    --background ${gene_list} \
    --deseq_results ${deseq_results}
  """
}
