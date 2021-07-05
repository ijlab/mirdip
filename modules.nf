process MIRDIP5_UPDATE_IDENTIFIERS {
  conda "/home/waddelld/.conda/envs/mirbaseconverter"
  publishDir "$params.updated_ids_folder"

  input:
    file(not_updated)

  output:
    file "${not_updated}.ids_updated.tsv"

  """
  Rscript $params.convert_mirbase_id_script -t $not_updated -n FALSE -o "${not_updated}.mirbase_v22_IDs.tsv" -s \$'\t' -c "2" > "${not_updated}.mirbase_VERSIONS" 2>&1
  awk -F \$'\t' '\$2 == "NA"' <"${not_updated}.mirbase_v22_IDs.tsv" | sort | uniq > "${not_updated}.mirbase_v22_IDs.tsv.dead.uniq"
  awk -F \$'\t' '\$1 != \$2' <"${not_updated}.mirbase_v22_IDs.tsv" | sort | uniq > "${not_updated}.mirbase_v22_IDs.tsv.changed.uniq"
  python $params.update_data_source_ids_script --path="${not_updated}" --checker "${params.hgnc_symbol_check_results}" --outpath="${not_updated}.ids_updated.tsv"
  """
}

process MIRDIP5_BENCHMARK_ONE_RESOURCE {
  conda "/home/tools/anaconda3/envs/mirbaseconverter"
  publishDir "$params.bench_collection"

  input:
    file(scores)

  output:
    file "${scores}.benchmarked.tsv"

  """
  Rscript $params.benchmark_script_path \
          -c `pwd` \
          -t $scores \
          -g "${params.gold_standard}" \
	  -s "${params.bench_style}"
  """
}

process MIRDIP5_INTEGRATE_AND_PLOT_BENCHMARKS {
  conda "/home/tools/anaconda3/envs/mirbaseconverter"

  input:
    val bench_collection

  output:
    file "${params.prec_recall_pdf}"

  """
  Rscript $params.integrate_and_plot_path \
          -c `pwd` \
          -d bench_collection \
	  -o "${params.prec_recall_pdf}" \
          -g "${params.gold_standard}"
  """
}

process MIRDIP5_RUN_NOISYOR {
  conda "/home/tools/anaconda3/envs/mirbaseconverter"

  input:
    file "*.benchmarked.tsv"

  output:
    file "${params.integrated_score_file}"

  """
  Rscript $params.noisyor_script \
          -c `pwd` \
          -d "${params.bench_collection}" \
          -o "${params.integrated_score_file}"
  """
}