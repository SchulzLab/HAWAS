nextflow.enable.dsl = 2

process EPIEXPRESS_PREDICTION {
  tag "EpiExpress_Prediction_${mode}"

  input:
  val epiexpress_repo_dir
  val bigwigs_input
  val cre_rf_model_folder_input
  val binned_cnn_model_folder_input
  val provided_input_dir
  val gene_file
  val mode
  val correlation_cutoff
  val cores_param

  output:
  path "epiexpress_predictions_output"

  script:
  def ABS_GENE_FILE = gene_file
  """
  #!/bin/bash
    set -e

  # unset local R library paths so conda env wins
  unset R_LIBS
  unset R_LIBS_USER

  # Setup
  WORK_DIR=\$(pwd)
  mkdir -p epiexpress_predictions_output
  
  cat <<EOF > "\${WORK_DIR}/epiexpress_input.json"
{
  "bigwigs": "${bigwigs_input}",
  "mode": "${mode}",
  "out_folder": "\${WORK_DIR}/epiexpress_predictions_output",
  "gene_file": "${ABS_GENE_FILE}",
  "cre_rf_model_folder": "${cre_rf_model_folder_input}",
  "binned_cnn_model_folder": "${binned_cnn_model_folder_input}",
  "provided_input": "${provided_input_dir}",
  "cores": ${cores_param},
  "correlation_cutoff": ${correlation_cutoff}
}
EOF

  # Clean previous CRE_input files to avoid gzip "already exists" errors
  if [ -d "\${WORK_DIR}/epiexpress_predictions_output/CRE_input" ]; then
      rm -rf "\${WORK_DIR}/epiexpress_predictions_output/CRE_input"/*
  fi

  echo "Step 1: BigWigToInput"
  python "${epiexpress_repo_dir}/Prediction/BigWigToInput.py" "\${WORK_DIR}/epiexpress_input.json"

  MODE="${mode}"
  echo "Step 2: predictions (mode: \$MODE)"
  
  if [[ "\$MODE" == "CRE-RF" ]]; then
      echo "Running CRE-RF prediction..."
      Rscript "${epiexpress_repo_dir}/Prediction/Prediction_CRE_RF_parallel.R" "\${WORK_DIR}/epiexpress_input.json"
      
  elif [[ "\$MODE" == "Binned-CNN" ]]; then
      echo "Running Binned-CNN prediction..."
      python3 "${epiexpress_repo_dir}/Prediction/Prediction_Binned_CNN.py" "\${WORK_DIR}/epiexpress_input.json"
      
  elif [[ "\$MODE" == "all" ]] || [[ "\$MODE" == "CRE-RF,Binned-CNN" ]]; then
      echo "Running both CRE-RF and Binned-CNN predictions..."
      Rscript "${epiexpress_repo_dir}/Prediction/Prediction_CRE_RF_parallel.R" "\${WORK_DIR}/epiexpress_input.json"
      python3 "${epiexpress_repo_dir}/Prediction/Prediction_Binned_CNN.py" "\${WORK_DIR}/epiexpress_input.json"
      
  else
      echo "ERROR: Unknown mode '\$MODE'"
      echo "Valid options: CRE-RF, Binned-CNN, all, CRE-RF,Binned-CNN"
      exit 1
  fi

  echo "EpiExpress prediction completed"
  """
}
