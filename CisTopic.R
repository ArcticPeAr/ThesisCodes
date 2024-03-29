#Load data from command line in snakemake
args <- commandArgs(trailingOnly = TRUE)


#To impute missing values in the data amelia is used
library("Amelia")
library("tidyverse")
library("readODS")

argCounter =1

NumberOfCores = args[argCounter]
NumberOfCores = as.integer(NumberOfCores)
argCounter = argCounter + 1
cat("\n","NumberOfCores: ", NumberOfCores, "\n")


#String for the BRCA - loop
whichLoopString <- args[argCounter]
argCounter = argCounter + 1

cat (" Sample: ", whichLoopString, "\n")
#Number of cores to be used

#The methylation table
CpG <- args[argCounter]
argCounter = argCounter + 1


#The coordinates of the CpGs with respect to the patients
coords <- args[argCounter]
argCounter = argCounter + 1
#metadata of the patients. Needs to be ordered
premeta <- args[argCounter]
argCounter = argCounter + 1


#methTable <- read.table(CpG, header=TRUE, sep="\t")
methTable <- readRDS(CpG)
Coords <- read.table(coords)

#Massage the Coords to be in the right format for methTable
mergedCol <- unite(Coords, comb, V2, V3, sep = "-", remove = TRUE, na.rm = FALSE)
mergedChrDone <- unite(mergedCol, Coords, V1, comb, sep = ":", remove = TRUE, na.rm = FALSE)
names(mergedChrDone)[names(mergedChrDone) == "V4"] <- "probes"
mergedChrDone$probes <- as.character(mergedChrDone$probes)

#Massage the methTable to be in the right format for Coords
names(methTable)[names(methTable) == "X"] <- "probes"
methTable$probes <- as.character(methTable$probes)

methTable <- left_join(mergedChrDone, methTable, by = c("probes"))

methTable <- within(methTable, rm("probes"))

methTable <- methTable[which(rowMeans(!is.na(methTable)) > 0.5), ]


# Follows is a loop for bounds to prevent AmeliaII to impute some negative values in the dataframe

nc = ncol(methTable)

loopyBound = rbind(c(2, 0, 1))

n=3
while (n <= nc)
{
loopyBound = rbind(loopyBound, c(n, 0, 1))
n=n+1
}

idvars = c('Coords')

AmCores = as.integer(floor(NumberOfCores/2))

#parallelization of Amelia needs more than 1 core. If only 1 core is available, then Amelia is run without parallelization
if (AmCores > 1) {
#multicore for Amelia
library(parallel)
cl <- makeCluster(AmCores)

options("Amelia.clustermode" = cl)

methTable <- amelia(methTable, m = 1, idvars = idvars,parallel = "multicore")

stopCluster(cl)
} else {
methTable <- amelia(methTable, m = 1, idvars = idvars)
}

methTable <- methTable$imputations[[1]]

methTable2 <- methTable[,-1]
rownames(methTable2) <- methTable[,1]


#The metadata of the patients is prepared
cat("Metadata is being prepared\n")
samplemeta <- read.table(premeta, header=TRUE, sep="\t")

if (whichLoopString == "BRCA-US")
{
metaframe <- samplemeta %>% select(donor_id, PAM50, ER.Status, PAM50_genefu, PR.Status, HER2.Final.Status)

#replace numbers as such: 1=Negative, 2=Positive, and LuminalA and LuminalB to "Luminal A" and "Luminal B"
metaframe$ER.Status <- gsub("1","Negative", metaframe$ER.Status)
metaframe$ER.Status <- gsub("2","Positive", metaframe$ER.Status)
metaframe$PAM50 <- gsub("LuminalA","Luminal A", metaframe$PAM50)
metaframe$PAM50 <- gsub("LuminalB","Luminal B", metaframe$PAM50)

#If redundancy in data (E.g same donors multiple times)
metaSD <- metaframe %>% distinct(donor_id, .keep_all = TRUE)

} else 
{ 
metaframe <- samplemeta %>% select(icgc_donor_id, Subtype_Selected, Subtype_Immune_Model_Based)
metaSD <- metaframe %>% distinct(icgc_donor_id, .keep_all = TRUE)
}  

##########Make sure the metadata is correct!

#set rownames
meta2 <- metaSD
rownames(meta2) <- metaSD[,1]


meta <- meta2 %>% select(-1)


###########RemoveNonMetaCpGRows############

MTList <- colnames(methTable2)
metaList <- rownames(meta)

for (item in MTList)
{
    if (item %in% metaList==TRUE)
    {
        MTList <- MTList[! MTList %in% c(item)]
    }
}

#MTList = list of rows needing to be added to metadata
rn <- row.names(meta)
meta[nrow(meta) + seq_along(MTList), ] <- NA 
row.names(meta) <- c(rn, MTList)


#replace NA with string
#meta <- meta %>% mutate_all(as.character)
meta <- meta %>% replace(is.na(.), "NA")

#Sort meta
meta <- meta[order(row.names(meta)),]

MTList <- colnames(methTable2)

#remove rows not in methTable2 (to be used in topicBina.py)
meta <- meta[row.names(meta) %in% MTList, ]


#free memory:
rm(methTable)
rm(Coords)
rm(samplemeta)
rm(metaframe)
rm(metaSD)


methTable2ods <- args[argCounter]
argCounter = argCounter + 1
#save the methTable2 as odt (because of ods has a more compact storage than csv)
write_ods(methTable2, "methTable2.ods")

######################################################################################
#CISTOPIC
######################################################################################
suppressWarnings(library("cisTopic"))

#Name the PDF after what cancer type and TF:
pdfpdf <- args[argCounter]
argCounter = argCounter + 1
PDFName <- pdfpdf


cisTopicObject <- createcisTopicObject(is.acc=0.5, min.cells=0, min.regions=0, count.matrix=data.frame(methTable2)) #Set is.acc=0 or is.acc=0.01 | min.cells=0, min.regions=0

cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = meta)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(4:18), seed=123, nCores=NumberOfCores, addModels=FALSE)

cat("cisTopicObject created\n")
pdf(PDFName, height=10, width=15)
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cat("Model selected\n")

ctoOut <- args[argCounter]
argCounter = argCounter + 1
saveRDS(cisTopicObject, file=ctoOut)

topicAssigToPatientOut <- args[argCounter]
argCounter = argCounter + 1
#### Write out topic assignments to the patients
write_ods(cisTopicObject@selected.model$document_expects, topicAssigToPatientOut)

RegScrPrtopicOut <- args[argCounter]
argCounter = argCounter + 1
write_ods(cisTopicObject@region.data, RegScrPrtopicOut)
#write.xlsx(cisTopicObject@region.data, file="/storage/mathelierarea/processed/petear/analysis/test/RegScrPrtopic.csv")

#### Unnormalized region assignments
RegAssigUnormalOut <- args[argCounter]
argCounter = argCounter + 1
write_ods(cisTopicObject@selected.model$topics, RegAssigUnormalOut)

#### Binarized region assignments
binaOut <- args[argCounter]
bina <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)
write_ods(bina, binaOut)
