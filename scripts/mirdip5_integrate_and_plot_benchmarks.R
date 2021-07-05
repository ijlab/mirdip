#!/usr/bin/env Rscript

# usage: mirdip5_integrate_and_plot_benchmarks.R -d /path/to/rda_files.Rda -o /output/precision_recall_plots.pdf

library(HGNChelper)
library(data.table)
library(miRBaseConverter)
library(docopt)
library(stringi)
library(utils)

'Pick up the individual statistics caluclated from each benchmarked resources and do some Precision/Recall curves, then print to names PDF file output.

Usage:
 mirdip5_integrate_and_plot_benchmarks.R -d=<DIRECTORY> -o=<OUT> -c=<WORKDIR> -g=<GOLD>

Options:
  -h --help             Show the help message.
  -c --wd=<WORKDIR>     Set the working directory to this.
  -d --dir=<DIRECTORY>  Path to the parent dir of the Rda files corresponding to stats from bench.
  -o --out=<OUT>        Name of output PDF with precision plots.
  -g --gs=<GOLD>        Gold standard data for benchmarking [default: /home/tokar/workspace/mirDIPv4/gold_standard/gold_standard_large.txt]
  
' -> doc

arguments <- docopt(doc, version = 'mirdip5_integrate_and_plot_benchmarks.R 0.1')
print(arguments)

setwd(arguments$wd)

# Read gold-standard
gs = fread(arguments$gs, sep = '\t', head = T)

# Take only benchmarking data
benchmark = subset(gs, Data.class == 'benchmark')[,1:2]
# setkey(benchmark, 'Gene symbol', 'miRNA')
# bench.miR = unique(benchmark$miRNA)
# bench.gen = unique(benchmark$`Gene symbol`)
setkey(benchmark, 'Gene symbol', 'miRNA')
bench.miR = unique(benchmark$miRNA)
bench.gen = unique(benchmark$`Gene symbol`)

calc.stats = function(r, p.data, p.bench){
  y = p.data[V8 < r & V8 > (r - dW) ,]
  m = median(y$V8)
  y = y[,1:2]
  n = nrow(unique(y))
  tp = nrow(fintersect(y, p.bench))
  fp = nrow(fsetdiff(y, p.bench))
  fn = nrow(fsetdiff(p.bench, y))
  precis = tp / (tp + fp)
  recall = tp / (tp + fn)
  fscore = 2 * ((precis * recall) / (precis + recall))

  s = c(r, m, n, tp, fp, fn, precis, recall, fscore)
  return(s)
}

load_stats <- function(rdata_file) {
  # fresh, empty namespace!
  load(rdata_file)
  my_vars <- ls()
  #stopifnot(length(my_vars) == 1)
  #stopifnot(my_vars == "stats")
  stats
}

resources = Sys.glob(file.path(arguments$dir, "*.bench.Rda"))
n = length(resources)
results = list()
dR = 0.001
dW = 0.01
Rs = seq(dW, 1., dR)
for (i in 1:n){

  # Read data
  print(resources[i])
  results[[i]] <- load_stats(resources[i])

}

save(results, file = './dr0001_dw_001_gliding_large_GS.Rdata')

rmse = function(x, y){
  d = cbind(x, y)
  d = d[complete.cases(d),]
  er = d[,1] - d[,2]
  sqrt(mean(er^2))
}

pdf(arguments$out)
for (i in 1:n){
  d = results[[i]]
  
  #d = d[which.max(d$Precision):nrow(d),]
  d = d[complete.cases(d),]
  #d = d[d$R > 0,]

  if (nrow(d) > 1){
    x = log(d$R)
    y = log(d$Precision)
    w = d$N
    m = lm(y ~ poly(x, 2))
    p = exp(predict(m, data.frame(x  = log(seq(0.01, 1, 0.01)))))

    plot(results[[i]]$R,
         results[[i]]$Precision,
         type = 'l',
	 xlab = 'Recall',
	 ylab = 'Precision',
         main = paste0("PR Curve for ", basename(resources[i])))

    lines(seq(0.01, 1, 0.01),
          p,
          col = 'red')

    #print(rmse(results[[i]]$Precision, p))
  }
}
dev.off()