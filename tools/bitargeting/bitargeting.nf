#!/usr/bin/env nextflow

mir_chunks = Channel.fromPath(params.mirs_file)
                    .splitFasta( by: params.mirs_per_process )
utr_chunks = Channel.fromPath(params.utrs_file)
	            .splitFasta( by: params.utrs_per_process )

process runBitargetingBatch {
  conda '/home/waddelld/rnatools/rnatools.yml'
  
  input:
  file 'chunk.fasta' from mir_chunks
  each file(utrs) from utr_chunks
  val app from params.bitargetingJar

  """
  python ${params.genconf_py} --bitargetingJar=${params.bitargetingJar} --micrornas=chunk.fasta --genes=$utrs --energy_filter=${params.energy_filter} --energy_cutoff=${params.energy_cutoff}
  java -jar $app config.txt parameters.txt
  """
}