#!/usr/bin/env nextflow

params.tgs_pl = "/home/waddelld/rnatools/targetscan/targetscan_70.pl"
params.mirs_file = "/home/waddelld/rnatools/all5_seqs.fa"
params.genes_file = "/home/waddelld/rnatools/aura/UTR_hg19.fasta"
params.mirs_per_process = 1

mir_chunks = Channel.fromPath(params.mirs_file)
                    .splitFasta( by: params.mirs_per_process )

// determine whether to use container for this one or not...
process runBitargetingBatch {
  conda '/home/waddelld/rnatools/rnatools.yml'
  cpus 1
  queue "dque"
  executor "slurm"
  memory "1 GB"
  time "10m"
  clusterOptions "--job-name=bitargeting_batch -N1"

  input:
  file 'chunk.fasta' from mir_chunks
  val app from params.tgs_pl
  val utrs from params.genes_file

  """
  perl $params.tgs_pl 
  """
}