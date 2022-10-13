# Post-processing of Contexts for RNAseq data for genes and miRNAs

Post-processing steps for mRNA as well as miRNA are available in the corresponding folders, `postprocessingMRNA` and `postprocessingMiRNA`. Each script in this folder belongs to one specific dataset, since the dataset is processed separately. All scripts are written in R and therefore you need to install R to run them. Moreover, they are named based on its dataset source, beginning with `dataCleanup` and ends with its source name  `dataCleanup<datasetname>.Rmd`. To ensure term consistency, we used the Disease Ontology and BRENDA Tissue Ontology to standardize context names. When a term was not present, its relationships were identified through other ontologies in the Ontology Lookup Service (OLS) - namely FMA, NCIT, UBERON, and OBA. Finally, we curated the remaining relationships to map them to terms already included. Contexts corresponding to cell lines and qualifiers outside normal and disease (for example, developmental stage) were not included in this release. Furthermore, all gene expression datasets were post-processed to ensure that all gene symbols were consistent with the HGNC-approved symbols.

# miRNA post-processing
The post-processing steps of miRNA expression values are described in the following:
1. Ensure that all miRNAs were updated to miRBase v.22 IDs (miRBase homo-sapiens data: `mirdip/mirdip-tissue/data/mirbase/mature_homo-sapiens_dataframe.txt`)
2. Quintize miRNA expression values to enable a more fine-grained analysis of miRNA abundance. Therefore, for each sample, any miRNA with an expression value of zero remained so. The remaining non-zero values were converted to a number between one and five that represented which of the 20th percentiles of non-zero values it corresponds to.
3. Average its quantile-normalized values per miRNA (biological replicates, for instance).
4. Ensure context names are standardized.
5. Convert all miRNA expression values into binary values.


# mRNA post-processing
The post-processing steps of mRNA expression values are described in the following:
1. Average mRNA expression values per gene (biological replicates, for instance).
2. Ensure context names are standardized.
3. Convert all expression values into binary values.
