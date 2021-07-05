#!/usr/bin/env Rscript

# usage: mirdip5_bench_one_file.R -t /path/to/table.tsv

library(HGNChelper)
library(data.table)
library(miRBaseConverter)
library(docopt)
library(stringi)

'Just bench a single input file so that we do not have to wait for all of them to run serially. Load the saved stats DF to do the next part when all finished.

Usage:
 mirdip5_bench_one_file.R -t=<TABLE>

Options:
  -h --help             Show the help message.
  -t --table=<TABLE>    Path to the TSV file were working with.
  -g --gold=<GOLD>      Path to the gold standard we are using. Expects Validation and Benchmark classes.
  
' -> doc