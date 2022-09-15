/*
 * Runs miranda and parses output into something we can live with.
 */
process BATCH_PREDICT_miranda {
  conda "/home/waddelld/.conda/envs/rnatools" 

  input:
    tuple file(mir), file(utrs)

  output:
    file 'miranda.hits'
    file 'miranda.scans'

  """
  $params.miranda_executable $mir $utrs -out miranda.out
  grep -E "^>[^>]" miranda.out | tr -d '>' > miranda.hits
  grep ">>" miranda.out | tr -d '>' > miranda.scans
  """
}
