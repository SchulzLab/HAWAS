nextflow.enable.dsl = 2

process GENERATE_MATRIX {
  tag "${output_filename}"
  publishDir "${output_dir_matrix}", mode: 'copy'

  input:
  path predictions_dir
  val output_dir_matrix
  val n_cores
  val output_filename

  output:
  path "${output_filename}", emit: matrix_file

  script:
  """
  Rscript ${projectDir}/bin/generate_matrix.R \
    --input_dir ${predictions_dir} \
    --output_file ${output_filename} \
    --n_cores ${n_cores}
  """
}
