#!/usr/bin/env Rscript

# usage: mirdip5_run_noisyOR.R -d /path/to/individual_scores -o mirdip5_integrated_scores.txt

library(HGNChelper)
library(data.table)
library(miRBaseConverter)
library(docopt)
library(stringi)
library(dplyr)
library(utils)

'Runs Noisy-OR integrated score calculation for preprocessed, standardized, and normalized resources of mirdip.

Usage:
 mirdip5_run_noisyOR.R -d=<DIRECTORY> -o=<OUT> -c=<WORKDIR>

Options:
  -h --help             Show the help message.
  -c --wd=<WORKDIR>     Set the working directory to this.
  -d --dir=<DIRECTORY>  Path to the parent dir of the individual datasets.
  -o --out=<OUT>        Name of output file with the integrated scores.
  
' -> doc

arguments <- docopt(doc, version = 'mirdip5_run_noisyOR.R 0.1')
print(arguments)

setwd(arguments$wd)

# note that the dir is expected to have only files formatted
# like the resources_final are.
resources = Sys.glob(file.path(arguments$dir, "*.benchmarked.tsv"))
noisy.or = function(x)1 - exp(sum(log(1 - x)))
n = length(resources)
results = list()
for (i in 1:n){
  results[[i]] = fread(resources[i], sep='\t', head=F)
}
results = rbindlist(results)
results = results[, N := .(.N), by = .(V1, V2)]
results = results[, Score := .(noisy.or(V4)), by = .(V1, V2)]
results = results[, Sources := .(list(V5)), by = .(V1, V2)]
setorder(results, -Score)
results = unique(results, by = c(1,2))
results = results[,c(1,2,8:10)]
setnames(results, c('gene', 'MicroRNA', 
                'N', 'Score', 'Sources'))
mirdip_vh = results %>% group_by(MicroRNA) %>% filter(Score >= quantile (Score, 0.99)) %>% mutate(Class_2 = "V")
mirdip_h = results %>% group_by(MicroRNA) %>% filter(Score >= quantile (Score, 0.95) & Score < quantile (Score, 0.99)) %>% mutate(Class_2 = "H")                      
mirdip_m = results %>% group_by(MicroRNA) %>% filter(Score >= quantile (Score, 0.66) & Score < quantile (Score, 0.95)) %>% mutate(Class_2 = "M")
mirdip_l = results %>% group_by(MicroRNA) %>% filter(Score < quantile (Score, 0.66)) %>% mutate(Class_2 = "L")
mirdip_all = rbind(mirdip_vh, mirdip_h, mirdip_m, mirdip_l)
# all scores with all classes
fwrite(mirdip_all,
       basename(arguments$out),
       sep = '\t',
       row.names = F,
       quote = F)
# very high
fwrite(mirdip_vh,
       paste0(basename(arguments$out), ".very_high_score_class.txt"),
       sep = '\t',
       row.names = F,
       quote = F)
# high
fwrite(mirdip_h,
       paste0(basename(arguments$out), ".high_score_class.txt"),
       sep = '\t',
       row.names = F,
       quote = F)
# med
fwrite(mirdip_m,
       paste0(basename(arguments$out), ".medium_score_class.txt"),
       sep = '\t',
       row.names = F,
       quote = F)
# low
fwrite(mirdip_l,
       paste0(basename(arguments$out), ".low_score_class.txt"),
       sep = '\t',
       row.names = F,
       quote = F)