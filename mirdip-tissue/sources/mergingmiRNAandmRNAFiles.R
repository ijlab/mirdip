###################################################################################
### Author: Anne-Christin Hauschild
### Date:   2022-08-29

###   Create Final Files   ####################################################
###################################################################################


###################################################################################
###   Libarys and Paths  ####################################################
###################################################################################
library(scales)
library(stringr)

root_path <- "/home/mi/ownCloud/mywork/projects/MirDIP5.0/"
setwd(paste(sep="", root_path, "src/"))

#Linux
source_path <<- "./"
root_path <<- "../"


## Define Standard folders
data_path <- paste(root_path, "data/", sep="")
results_path <- paste(root_path, "results/", sep="")

##Inputfiles:
filesmiRNA<- paste(data_path,"miRNA_final",sep="")
filesmRNA<- paste(data_path,"RNA_final",sep="")
fileMiDIPgenes<- paste(data_path, "mirDIP_available_genes_U.txt", sep="")
filehgncSymbolCheck<- paste(data_path, "hgnc-symbol-check.csv", sep="")



## Outputfiles
allTissuesMirsFile<-paste(results_path, "allTissues_mirs_220904.txt", sep="")
allMirsFile<-paste(results_path, "all_mirs_220904.txt", sep="")
allTissuesGenesFile <- paste(results_path, "allTissues_genes_220904.txt", sep="")
allGenesFile<-paste(results_path, "all_genes_220904.txt", sep="")
allSourcesFile<-paste(results_path, "all_sources_220904.txt", sep="")
allmirnaExpressionFile<-paste(results_path, "all_miRNA_Expression_220904.txt", sep="")
allmrnaExpressionFile<-paste(results_path, "all_mRNA_Expression_220904.txt", sep="")


#########################################################################
############ Read Check Files:
#########################################################################
validGenes<- read.csv(fileMiDIPgenes, sep="", header = TRUE)
dim(validGenes)
head(validGenes)

geneCheck<- read.csv(filehgncSymbolCheck, sep=",", header = TRUE)
dim(geneCheck)
head(geneCheck)
unique(geneCheck[,"Match.type"])
#geneCheck<-geneCheck[geneCheck[,"Match.type"]!="Unmatched" ,]
geneMapping<- geneCheck[,"Approved.symbol"]
names(geneMapping)<- geneCheck[,"Input"]
head(geneMapping)

#########################################################################
############ Process miRNA Tissues Files:
#########################################################################

#Get list of miRNA data files
longdataFilesMIR <- list.files(filesmiRNA, pattern="*longData.txt")
longdataFilesMIR
#

# Init variables
allTissuesMirs<-c()
allMirs<-c()
allSources<-c()
longDataMerged<-c()

#Iterate and load all miRNA data files
longfile<-longdataFilesMIR[12]
for(longfile in longdataFilesMIR){
  
  print(longfile)
  #Read data
  datamiRNA<- read.csv(paste(filesmiRNA,"/", longfile, sep="") , sep="", header = TRUE)
  colnames(datamiRNA) <- c("miR", "Seq","Canonical","Source", "Tissue","Scale", "Binary")
  head(datamiRNA)
  # Get source, miRNAs, tissues
  uSource<-unique(datamiRNA[,"Source"] )
  uMirs<-unique(datamiRNA[,"miR"] )
  allMirs<-c(allMirs, uMirs)
  uTissuesMirs<-unique(datamiRNA[,"Tissue"] )
  allTissuesMirs<-c(allTissuesMirs, uTissuesMirs)
  
  
  allSources<-c(allSources, uSource)
  longDataMerged<-rbind(longDataMerged, datamiRNA)
  print(paste(uSource, "#Mirs", length(uMirs), "#Tissues",  length(uTissuesMirs)))
#  print(head(uTissuesMirs))
  
  #reset temp variables
  uSource<-c()
  uMirs<-c()
  uTissuesMirs<-c()
  datamiRNA<-c()
}
  
#Write required miRNA files for MirDIP
allTissuesMirs<-matrix(data = unique(allTissuesMirs), ncol = 1)
write.table(allTissuesMirs,  file=allTissuesMirsFile, sep="\t", row.names=FALSE, quote=FALSE)

allMirs<-matrix(data = unique(allMirs), ncol = 1)
write.table(allMirs,  file=allMirsFile, sep="\t", row.names=FALSE, quote=FALSE)

allSources<-matrix(data = c(allSources, rep("Yes", length(allSources))), ncol = 2, nrow=length(allSources) )
colnames(allSources)<-c("Dataset", "Scale")
write.table(allSources,  file=allSourcesFile, sep="\t", row.names=FALSE, quote=FALSE)

write.table(longDataMerged,  file=allmirnaExpressionFile, sep="\t", row.names=FALSE, quote=FALSE)




#########################################################################
############ Process mRNA Tissues Files:
#########################################################################

#Get list of mRNA data files
longdataFilesMRNA <- list.files(filesmRNA, pattern="*longData.txt")
longdataFilesMRNA


# Init variables
allTissuesGenes<-c()
allGenes<-c()
longDataMerged<-c()

#Iterate and load all mRNA data files
#longfile<-longdataFilesMRNA[7]
for(longfile in longdataFilesMRNA){
  #Read data
  datamRNA<- read.csv(paste(filesmRNA,"/", longfile, sep="") , sep="", header = TRUE)
  head(datamRNA)
  #Standardize columns
  colnames(datamRNA) <- c("Gene", "Source", "Tissue", "Binary")
  dim(datamRNA)
  head(datamRNA)
  #Use "hgnc-symbol-check" to select known Gene Symbols
  mappedGenes <- geneMapping[datamRNA[,"Gene"]]
  datamRNA[,"Gene"] =  mappedGenes
  #Remove unmatched genes
  idx<- unlist(mappedGenes!="" & !is.na(mappedGenes!=""))
  datamRNA<- datamRNA[idx,]

  # Get source, genes, tissues
  uSource<-unique(datamRNA[,"Source"] )
  uGenes<-sort(unique(datamRNA[,"Gene"] ))
  uTissuesGenes<-unique(datamRNA[,"Tissue"] )
  allGenes<-c(allGenes, uGenes)
  allTissuesGenes<-c(allTissuesGenes, uTissuesGenes)

  print(paste(uSource, "#Mirs", length(uGenes), "#Tissues",  length(uTissuesGenes)))
  print(dim(datamRNA))
  print(head(datamRNA))
  longDataMerged<-rbind(longDataMerged, datamRNA)

  #reset temp variables
  uSource<-c()
  uGenes<-c()
  uTissuesGenes<-c()
  datamRNA<-c()
}

#Write required mRNA files for MirDIP
allTissuesGenes<-matrix(data = unique(allTissuesGenes), ncol = 1)
write.table(allTissuesGenes,  file=allTissuesGenesFile, sep="\t", row.names=FALSE, quote=FALSE)

allGenes<-matrix(data = sort(unique(allGenes)), ncol = 1)
write.table(allGenes,  file=allGenesFile, sep="\t", row.names=FALSE, quote=FALSE)

write.table(longDataMerged,  file=allmrnaExpressionFile, sep="\t", row.names=FALSE, quote=FALSE)

