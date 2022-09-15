# mirdip5 pipeline

Runs according to `nextflow.config` configuration file, which controls things such as the source data directory, the benchmark/evaluation pairs, etc. It is crucial to update the `params` block of this configuration file in full to ensure correct results:

```groovy
params {
  scores_path = "/home/waddelld/rnatools/mirdip5/resources_map_these_ids"
  convert_mirbase_id_script = "/home/waddelld/rnatools/mirdip5/scripts/convert_mirbase_ids.R"
  update_data_source_ids_script = "/home/waddelld/rnatools/mirdip5/scripts/mirdip5_update_data_source_ids.py"
  hgnc_symbol_check_results = "/home/waddelld/rnatools/mirdip5/old_and_mirnatip_mirzag_targetscan_id_mapping.csv"
  updated_ids_folder = "/home/waddelld/rnatools/mirdip5/ids_updated"
  update_gene_symbols_script = "/home/waddelld/rnatools/mirdip5/scripts/update_gene_symbols.py"
  benchmark_script_path = "/home/waddelld/rnatools/mirdip5/scripts/mirdip5_bench_one_file.R"
  bench_collection = "/home/waddelld/rnatools/mirdip5/benchmarks_OLD_AND_ID_ISSUE"
  bench_style = "benchmark"
  gold_standard = '/home/waddelld/rnatools/mirdip5/plat_large_three_cols_only.tsv'
  integrate_and_plot_script = "/home/waddelld/rnatools/mirdip5/scripts/mirdip5_integrate_and_plot_benchmarks.R"
  prec_recall_pdf = '/home/waddelld/rnatools/mirdip5/mirdip5_benchmarks.pdf'
  noisyor_script = "/home/waddelld/rnatools/mirdip5/scripts/mirdip5_run_noisyOR.R"
  integrated_score_file = "mirdip5_integrated_scores.txt"
}
```

## How to run the pipeline

The following command will run the pipeline on `ijcluster`, produce an HTML report and timeline, as well as a `trace.txt` file with similar information about resource usage:

```bash
nextflow run map_ids.nf -profile ijcluster -with-report -with-timeline -with-trace
```

The ID mapping process relies on downloaded files `hgnc_complete_set.txt` as well as downloads from Ensembl to function correctly. Scripts in the `scripts/` directory will all depend on the libraries represented in the `mirbaseconverter.yml` conda environment file. Once that environment has been created and activated, you may inspect the options of each script and run each script independently simply by running it on the command line with the `-h` or `--help` options.

# running Noisy-OR integration script

The following command run on the results that are placed in the `params.publishDir` directory after successful execution of the above pipeline, will produce the integrated score. The `-d` option below is the path to `params.publishDir`, and the `-o` option is the path of the integrated score file.

```bash
# you may construct this environment from the mirbaseconverter.yml file in this repository
conda activate mirbaseconverter;
Rscript scripts/mirdip5_run_noisyOR.R \
        -c `pwd` \
        -d ./benchmarks_platinum_large_nodups/ \
        -o mirdip5_noisyor_final.txt
```
