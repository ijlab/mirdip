#---
#  title: "Processing Rahman expression data"
#  author: "Zuhaib Ahmed, Gitta Ekaputeri, Anne-Christin Hauschild"
#  date: "09/08/2021"
#  update: "29/08/2022" 
#  output: html_document
#---
  
setwd("/data/miRNA/schulzeData/")

fls <- list.files(pattern = "GSM")
head(fls)
data <- lapply(fls, function(x) {
  z <- read.table(x, header = T, sep = "\t")
  z <- z[,2:3]
  names(z) <- c("miR", str_split(x, "_")[[1]][1])
  return(z[1:2566,]) # The lines after 2566 are all empty for some reason.
})
names(data) <- sapply(str_split(fls, "_"), "[", 1)

head(data[[1]])
# data2 <- lapply(data, function(x) return(x[,2]))
# data2 <- as.data.frame(do.call(cbind, data2))
# data2 <- cbind(miR = data[[1]]$miR, data2)
data2 <- Reduce(function(x,y) merge(x, y, by = 1, all = T), data)
data2 <- data2[-1,] # For some reason, the first row was all NAs except for the GSM3015096 column. I'm deleting the first row
data2[1:5,1:5]

updte <- miRNAVersionConvert(data2$miR)
head(updte)
dim(updte)

mrbse <- read.table("../mirbase/mature_homo-sapiens_dataframe.txt", header = T, sep = "\t")
#mrbse <- read.table("../mirbase/mature_homo-sapiens_dataframe.txt", header = T, sep = "\t")
mrbse <- mrbse[-c(1,2),] # Remove the first two rows that don't match up with mirbase and are the only duplicates, and are the only sequences with T's for some reason.
head(mrbse)

mrbse2 <- merge(mrbse, updte, by.x = 1, by.y = 2, all.x = F, all.y = T)
head(mrbse2)
length(which(is.na(mrbse2$Name)))
mrbse2 <- mrbse2[-which(is.na(mrbse2$Name)),]
mrbse2 <- mrbse2[,-4]
head(mrbse2)
dim(mrbse2)



data3 <- merge(mrbse2, data2, by.x = 3, by.y = 1, all.x = T, all.y = F)
data3 <- data3[,-1] # Remove the original (old) IDs
data3[1:5,1:5]

dim(data3)
