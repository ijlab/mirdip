/*
 * Runs RNAhybrid.
 */
process BATCH_PREDICT_RNAhybrid {
  conda "/home/tools/anaconda3" 

  input:
    tuple file(mirs), file(utrs)

  output:
    file "${mirs.baseName}_${utrs.baseName}.out"

  """
  $params.rnahybrid_exe -c -s 3utr_human -t $utrs -q $mirs > "${mirs.baseName}_${utrs.baseName}.out"
  """
}
