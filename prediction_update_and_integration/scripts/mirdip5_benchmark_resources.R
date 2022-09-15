#!/usr/bin/env Rscript

# usage: mirdip5_bench_one_file.R -t /path/to/table.tsv

library(HGNChelper)
library(data.table)
library(miRBaseConverter)
library(docopt)
library(stringi)

'Just bench a single input file so that we do not have to wait for all of them to run serially. Load the saved stats DF to do the next part when all finished.

Usage:
 mirdip5_bench_one_file.R -d=<DIRECTORY>

Options:
  -h --help             Show the help message.
  -d --dir=<DIRECTORY>    Path to the parent dir of the Rda files corresponding to stats from bench.
  
' -> doc

arguments <- docopt(doc, version = 'mirdip5_bench_one_file.R 0.1')
print(arguments)

setwd('/gpfs/lb/mirdip5/tools/noisy_or')

# Read gold-standard
gs = fread('/home/tokar/workspace/mirDIPv4/gold_standard/gold_standard_large.txt', sep = '\t', head = T)

# Take only benchmarking data
benchmark = subset(gs, Data.class == 'benchmark')[,1:2]
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

pth = '/gpfs/lb/mirdip5/tools/noisy_or/mirdip5_resources_normalized/'
resources = list.files(pth, full.names = F)
n = length(resources)
results = list()
dR = 0.001
dW = 0.01
Rs = seq(dW, 1., dR)
for (i in 1:n){

  # Read data
  print(resources[i])
  fname = paste0(pth, resources[i])
  d = fread(fname, sep = '\t', head = F)
  d[,V8 := rank(d$V4)/nrow(d)]
  
  # Take subset overlaping with GS
  x = d[V1 %in% bench.gen | V2 %in% bench.miR,]
  colnames(x)[1:2] = colnames(benchmark)
  
  # Calculate statistics
  stats = lapply(Rs, calc.stats, x, benchmark)
  stats = do.call(rbind.data.frame, stats)
  colnames(stats) = c('r', 'R', 'N', 'TP', 'FP', 'FN', 'Precision', 'Recall', 'F-score')
  stats = na.omit(stats)
  
  
  x = log(stats$R)
  y = log(stats$Precision)
  m = lm(y ~ poly(x, 2))
  stats$Precision.fit = exp(predict(m, data.frame(x  = log(stats$R))))

  # EDIT: Don't do this, there is no reason to loop this serially.
  #       Save to a file, then for i do results[[i]]=stats and whatever
  #       else from there.
  # Add to results
  # results[[i]] = stats
  
  # Add benchmarked score
  d$V9 = exp(predict(m, data.frame(x = log(d$V8))))
  d$V9[d$V9 > max(stats$Precision)] = max(stats$Precision)
  
  # Reduce to unique miRNA-gene pairs
  setorder(d, -V9)
  d = unique(d, by = c(1,2))
  
  # Reorder columns
  d = d[,c(1, 2, 8, 9, 5:7)]
  
  # Write down
  fname = paste0('./mirdip5_resources_benchmarked/', resources[i])
  fwrite(d, fname, sep = '\t', col.names = F, row.names = F, quote = F)
}

save(results, file = './dr0001_dw_001_gliding_large_GS.Rdata')

rmse = function(x, y){
  d = cbind(x, y)
  d = d[complete.cases(d),]
  er = d[,1] - d[,2]
  sqrt(mean(er^2))
}

#load(file = './dr0005_dw_02_gliding_regular_GS.Rdata')
pdf('./precision_graphs_dr0001_dw001_gliding_large_DS.pdf')
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
         type = 'b',
         main = resources[i])

    lines(seq(0.01, 1, 0.01),
          p,
          col = 'red')

    #print(rmse(results[[i]]$Precision, p))
  }
}
dev.off()