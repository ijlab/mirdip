# Post-processing of MirDIP Tissue Contexts for Genes and miRNAs

The raw datasets collected were submitted to the Nextflow pipelines for further processing, using the Nextflow pipeline RNA-seq v. 3.8.1 (https://nf-co.re/rnaseq/3.8.1) for RNA seq data and smRNA-seq v. 2.0.0 (https://nf-co.re/smrnaseq) for miRNA set data. The datasets are then post-processed separately to ensure term consistency across the datasets. More detailed post-processing steps can be found in `mirdip/mirdip-tissue/sources/`. The `data/mirbase` folder containts the miRBase release used for the processing.
