nextflow.enable.dsl = 2

include {
  MIRDIP5_BENCHMARK_ONE_RESOURCE;
  MIRDIP5_INTEGRATE_AND_PLOT_BENCHMARKS;
  MIRDIP5_RUN_NOISYOR } from './modules.nf'
  
workflow bench_mirdip5_resources {
  take: mir_gene_scores_norm
  main:
    MIRDIP5_BENCHMARK_ONE_RESOURCE(mir_gene_scores_norm)
}

workflow integrate_and_plot_benchmarks {
  take: benchmarked_files
  main:
    MIRDIP5_INTEGRATE_AND_PLOT_BENCHMARKS(benchmarked_files)
}

workflow run_noisyor {
  take: benchmarked_files
  main:
    MIRDIP5_RUN_NOISYOR(benchmarked_files)
}

workflow {
    bench_mirdip5_resources(channel.fromPath("${params.scores_path}/*"))
}