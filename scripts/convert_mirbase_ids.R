#!/usr/bin/env Rscript

# usage: convert_mirbase_ids.R -t /path/to/table.tsv -c column_name_or_num_with_mirs -o output_table_name

library(miRBaseConverter)
library(docopt)
library(stringi)

'Requires miRBaseConverter and docopt packages installed. Something like conda activate mirbaseconverter could work. Will attempt to map the specified column in the table to mirbase version 22 by default, modify the -v parameter to choose a different target version. A summary of the best matches to all mirbase versions given the input list of mirbase IDs is printed by default - disable this by setting the -m paramter to FALSE if desired. The output of this script is a simple lookup table with unique mirbase IDs that were found in the input column from the input table, if any, in the first column - and in the second and last column, the corresponding mapping to the target mirbase version.

Usage:
  convert_mirbase_ids.R -t=<TABLE> -o=<OUTFILE> -s=<SEP> -c=<COL> [-n=<HEADER>] [-v=<TARGET>]

Options:
  -h --help             Show the help message.
  -t --table=<TABLE>    Path to the TSV file were working with.
  -s --sep=<SEP>        Separator that will be used in table files.
  -n --header=<HEADER>  Whether first line describes cols [default: FALSE]
  -c --col=<COL>        Column name or position for mirbase IDs
  -o --out=<OUTFILE>    Output TSV file, simply a map of distinct mirbase vX -> mirbase vY
  -v --target=<TARGET>  Target mirbase version to map to [default: 22]
  
' -> doc

arguments <- docopt(doc, version = 'convert_mirbase_ids.R 0.1')
print(arguments)

# read in the table.tsv file specified
matrix <- read.delim(arguments$table, sep = arguments$sep, header = ifelse(arguments$header %in% c("0", "TRUE"), TRUE, FALSE))
#print(sprintf("read matrix as length %1$s, with columns:\n%2$s", nrow(matrix), colnames(matrix)))
# if it's a digit, use it as the column index. Otherwise treat as name.
#print(arguments$col)
#print(str(arguments$col))
colname = ifelse(stri_detect_regex(arguments$col, "^[[:digit:]]+$"), paste("V", arguments$col, sep=''), arguments$col)
#print(colname)
mirs=matrix[[colname]]
#print(stri_detect_regex(arguments$col, "^[[:digit:]]+"))
print(sprintf("micrornas: %1$s", mirs))
# 4.2 says there could be more exact matches this way, redirect stdout/err here if you want the % for each version
#print(mirs)
version=checkMiRNAVersion(mirs, verbose=TRUE)
#print(paste("MicroRNAs matched mirbase version: ", version))
# this is the best matched version
# result1=miRNA_NameToAccession(mirs, version=version)
# grab the names from the accessions - seems it is the logic in the above/below function that adds mappings? Checked source code, meh.
# result2=miRNA_AccessionToName(result1[,2],targetVersion=paste0("v", arguments$target))

# yeah... that didn't work. We'll come back to it. For now, this is an updated mapping we can use.
converted_mirs=miRNAVersionConvert(mirs, targetVersion=paste0("v", arguments$target), exact=TRUE)
write.table(converted_mirs, file = arguments$out, sep = '\t', row.names = F, quote = F)