/*
 * Runs tarpmir.
 */
process BATCH_PREDICT_TarPmiR {
  conda "/home/waddelld/.conda/envs/tarpmir" 

  input:
    tuple file(mirs), file(utrs)

  output:
    file "${mirs}_${utrs}.bp"

  """
  python $params.prediction_script_path -a $mirs -b $utrs -m $params.prediction_model_path -p $params.prob
  """
}
