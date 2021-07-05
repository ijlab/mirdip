/*
 * predict_multi version uses suboptimal strategy for parallelization.
 * Instead, we use the signle miR, multiple UTR version and do our own
 * job scheduling.
 */
process BATCH_PREDICT_PITA {
  conda "/home/tools/anaconda3/"

  input:
    tuple file(mir), file(utrs)

  output:
    file '1_pita_results.tab'
    file '1_pita_results_targets.tab'

  """
  perl $params.prediction_script_path \
    -utr $utrs \
    -mir $mir \
    -prefix "${params.prefix}"
  """
}
