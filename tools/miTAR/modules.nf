/*
 * predict_multi version uses suboptimal strategy for parallelization.
 * Instead, we use the signle miR, multiple UTR version and do our own
 * job scheduling.
 */
process BATCH_PREDICT_miTAR {
  conda "/home/waddelld/.conda/envs/mitar" 

  input:
    tuple file(mir), file(utrs)

  output:
    file 'result.fa' 

  """
  ln -s "${params.mitar_utils_for_import}" utils.py
  python $params.prediction_script_path \
    -i1 $utrs \
    -i2 $mir \
    -o result.fa \
    -p $params.prob \
    -ns $params.numsites \
    -s $params.step \
    -m "${params.prediction_model_path}"
  """
}
