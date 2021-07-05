#!/usr/bin/env nextflow

params.mirmap_py = "/home/waddelld/rnatools/mirmap/miRmap-1.1/scripts/mirmap_runs.py"
params.mirs_file = "/gpfs/lb/mirdip5/mirbase_human_plus_78_from_mirgenedb.fa"
params.utrs_file = "/gpfs/lb/mirdip5/utrs/utrs_hg38_removed_unavailables_and_double_dash.fa"
params.mirs_per_process = 1
params.utrs_per_process = 100000

mir_chunks = Channel.fromPath(params.mirs_file)
                    .splitFasta( by: params.mirs_per_process )		    
utr_chunks = Channel.fromPath(params.utrs_file)
	            .splitFasta( by: params.utrs_per_process )

process runMirMapBatch {
  conda '/home/waddelld/.conda/envs/rnatools'
  cpus 1
  queue "dque"
  executor "slurm"
  memory "512 MB"

  input:
  file 'chunk.fasta' from mir_chunks
  each file(utrs) from utr_chunks
  val app from params.mirmap_py

  """
  export PYTHONPATH=/home/waddelld/rnatools/mirmap/miRmap-1.1/src
  python $app -a chunk.fasta -f $utrs -o mirmap_output.out -l /home/waddelld/rnatools/mirmap/miRmap-1.1/libs/lib-archlinux-x86_64
  """
}