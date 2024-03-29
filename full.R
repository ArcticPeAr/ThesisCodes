#Number of cores to be used
NoCores = 4
#Creates proper matrix before using topic modeling and performing great analysis
library("Amelia")

args <- commandArgs(trailingOnly = TRUE)

CpG <- args[1]
cat("CpG: ", CpG, "")
coords <- args[2]
premeta <- args[3]

methTable <- readRDS(CpG)
#methTable <- head(methTable, 10000)
Coords <- read.table(coords)
library(tidyverse)

#methTable <- read.table(snakemake@input[["CpG"]], header=TRUE, sep="\t")
#Coords <- read.table(snakemake@input[["coords"]])
#methTable <- read.table("/storage/mathelierarea/processed/petear/SnakemakeInputFiles/BRCA-US_methylation450.sorted.perDonor.tsv", header=TRUE, sep="\t")
#Coords <- read.table("/storage/mathelierarea/processed/petear/SnakemakeInputFiles/BRCA-US_FULL.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed")

mergedCol <- unite(Coords, comb, V2, V3, sep = "-", remove = TRUE, na.rm = FALSE)
mergedChrDone <- unite(mergedCol, Coords, V1, comb, sep = ":", remove = TRUE, na.rm = FALSE)
names(methTable)[names(methTable) == "X"] <- "probes"
methTable$probes <- as.character(methTable$probes)


names(mergedChrDone)[names(mergedChrDone) == "V4"] <- "probes"
mergedChrDone$probes <-as.character(mergedChrDone$probes)

methTable <- left_join(mergedChrDone, methTable, by = c("probes"))

methTable <- within(methTable, rm("probes"))

methTable <- methTable[which(rowMeans(!is.na(methTable)) > 0.5), ]


### Follows is a loop for bounds to prevent AmeliaII to make some negative values in the dataframe
nc = ncol(methTable)

loopyBound = rbind(c(2, 0, 1))

n=3
while (n <= nc)
{
loopyBound = rbind(loopyBound, c(n, 0, 1))
n=n+1
}

idvars = c('Coords')

methTable <- amelia(methTable, m = 1, idvars = idvars)

methTable <- methTable$imputations[[1]]

methTable2 <- methTable[,-1]
rownames(methTable2) <- methTable[,1]


#############PrepMeta###############

samplemeta <- read.table(premeta, header=TRUE, sep="\t")
#samplemeta <- read.table("/storage/mathelierarea/processed/petear/SnakemakeInputFiles/Meta/sampleinfo_TCGA_RNA_seq_cluster.txt", header=TRUE, sep="\t")

###OTHERS THAN BRCA::: _Subtype_mRNA', 'Subtype_Selected' & 'Subtype_Immune_Model_Based
#metaframe <- samplemeta %>% select(icgc_donor_id, Subtype_Selected, Subtype_Immune_Model_Based)


###BRCA_FULL
metaframe <- samplemeta %>% select(donor_id, PAM50, ER.Status, PAM50_genefu, PR.Status, HER2.Final.Status)



#replace numbers as such: 1=Negative, 2=Positive, and LuminalA and LuminalB to "Luminal A" and "Luminal B"
metaframe$ER.Status <- gsub("1","Negative", metaframe$ER.Status)
metaframe$ER.Status <- gsub("2","Positive", metaframe$ER.Status)
metaframe$PAM50 <- gsub("LuminalA","Luminal A", metaframe$PAM50)
metaframe$PAM50 <- gsub("LuminalB","Luminal B", metaframe$PAM50)

#If redundancy in data (E.g same donors multiple times)
#metaSD <- metaframe %>% distinct(icgc_donor_id, .keep_all = TRUE)
metaSD <- metaframe %>% distinct(donor_id, .keep_all = TRUE)

#set rownames
meta2 <- metaSD
rownames(meta2) <- metaSD[,1]


meta <- meta2 %>% select(-1)


###########RemoveNonMetaCpGRows############
#Cange to add metadata rows with NA for CpGs without metadata
#Find CpGs without row in meta
#
#som e i col men ikke i row.

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

###########################################
############CISTOPIC########
###########################################
suppressWarnings(library("cisTopic"))

#Name the PDF after what cancer type and TF:
pdfpdf <- args[4]
PDF.name <- pdfpdf

methTable2GZ <- args[5]
write.csv(methTable2, file=methTable2GZ)

#PDF.name <- snakemake@output[[1]]
#write.csv(methTable2, file=gzfile(snakemake@output[[2]]))
#write.csv(methTable2, file=gzfile("/storage/mathelierarea/processed/petear/analysis/test/methtable.csv.xz"))

#If removing NAs by deletion: (CHECK IF SHOULD USE AMPUTATION!!!
#cpg <- CpG.sample.tab[complete.cases(CpG.sample.tab), ]

cisTopicObject <- createcisTopicObject(is.acc=0.5, min.cells=0, min.regions=0, count.matrix=data.frame(methTable2)) #Set is.acc=0 or is.acc=0.01 | min.cells=0, min.regions=0

cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = meta)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(4:18), seed=123, nCores=NoCores, addModels=FALSE)
#cisTopicObject <- runModels(cisTopicObject, topic=c(8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60), seed=123, nCores=20, addModels=FALSE)
#cisTopicObject <- runModels(cisTopicObject, topic=c(2,4,6,8), seed=123, nCores=20, addModels=FALSE)

pdf(PDF.name, height=10, width=15)
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
#¤cisTopicObject <- selectModel(cisTopicObject)


##########################################################################################
#Save cisTopicObject
##########################################################################################
ctoOut <- args[69]
saveRDS(cisTopicObject, file=ctoOut)
#saveRDS(cisTopicObject, file="/storage/mathelierarea/processed/petear/analysis/test/cto.rds")

########## fICA ##############################################################################################################################################################################################################################
library(fastICA)

t_cisTo <- t(cisTopicObject@selected.model$document_expects)
donorMat  <- matrix(rownames(t_cisTo))
fica <- fastICA(t_cisTo, n.comp= 2, row.norm=FALSE, method = "C")


ficaS <- fica$S

colnames(ficaS) <- c("V1", "V2")
rownames(ficaS) <- c(donorMat)

metaficaS <- meta[rownames(meta) %in% rownames(ficaS), ]

metaficaS <- as.data.frame(metaficaS)
ficaS <- as.data.frame(ficaS)
ficaSMerge <- merge(metaficaS, ficaS, by = 0)

library(RColorBrewer)
################# GGPLOT2 for PAM50 ################
ficaSMergePAMSansNA <- ficaSMerge[!(ficaSMerge$PAM50=="NA"),]

figficaSPAM <- ggplot(ficaSMerge, aes(V1, V2))
figficaSPAM + geom_point(aes(colour = factor(PAM50)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSPAMSansNA <- ggplot(ficaSMergePAMSansNA, aes(V1, V2))
figficaSPAMSansNA + geom_point(aes(colour = factor(PAM50)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

################# GGPLOT2 for ER.Status ################

ficaSMergeERSansNA <- ficaSMerge[!(ficaSMerge$ER.Status=="NA"),]

figficaSER <- ggplot(ficaSMerge, aes(V1, V2))
figficaSER + geom_point(aes(colour = factor(ER.Status)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSERSansNA <- ggplot(ficaSMergeERSansNA, aes(V1, V2))
figficaSERSansNA + geom_point(aes(colour = factor(ER.Status)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())


################# GGPLOT2 for PAM50_genefu ################
ficaSMergePGSansNA <- ficaSMerge[!(ficaSMerge$PAM50_genefu=="NA"),]

figficaSPG <- ggplot(ficaSMerge, aes(V1, V2))
figficaSPG + geom_point(aes(colour = factor(PAM50_genefu)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSPGSansNA <- ggplot(ficaSMergePGSansNA, aes(V1, V2))
figficaSPGSansNA + geom_point(aes(colour = factor(PAM50_genefu)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())


################# GGPLOT2 for PR.Status ################
ficaSMergePRSansNA <- ficaSMerge[!(ficaSMerge$PR.Status=="NA"),]

figficaSPR <- ggplot(ficaSMerge, aes(V1, V2))
figficaSPR + geom_point(aes(colour = factor(PR.Status)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSPRSansNA <- ggplot(ficaSMergePRSansNA, aes(V1, V2))
figficaSPRSansNA + geom_point(aes(colour = factor(PR.Status)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())


################# GGPLOT2 for HER2.Final.Status ################
ficaSMergeHER2SansNA <- ficaSMerge[!(ficaSMerge$HER2.Final.Status=="NA"),]

figficaSHER2 <- ggplot(ficaSMerge, aes(V1, V2))
figficaSHER2 + geom_point(aes(colour = factor(HER2.Final.Status)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSHER2SansNA <- ggplot(ficaSMergeHER2SansNA, aes(V1, V2))
figficaSHER2SansNA + geom_point(aes(colour = factor(HER2.Final.Status)), size = 5) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())


######################################################
#####Umap and further cisTopic
######################################################
###Remove warnings
oldw <- getOption("warn")
options(warn = -1)
####Warning message:
#failed creating initial embedding; using random embedding instead




cisTopicObject <- runUmap(cisTopicObject, target ='cell')
cellTopicHeatmap(cisTopicObject, col.low="blue", col.mid="white",col.high="red", colorBy=c("PAM50","ER.Status","PAM50_genefu","PR.Status","HER2.Final.Status"))
cellTopicHeatmap(cisTopicObject, method = "Probability",col.low="blue", col.mid="white",col.high="red", colorBy=c("PAM50","ER.Status","PAM50_genefu","PR.Status","HER2.Final.Status"))



plot.new()
mtext("regionScore by NormTop")


cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
#par(mfrow=c(3,5))
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)


cisTopicObject <- GREAT(cisTopicObject, genome='hg19', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
ontologyDotPlot(cisTopicObject, top=5, var.y='name', order.by='Binom_Adjp_BH')


# plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
#      xaxt='n', yaxt='n', xlab='', ylab='')
# topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
#text(1,4, getwd() , pos=4)

cisTopicObject <- runUmap(cisTopicObject, target="cell", n_components=3)
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("PAM50","ER.Status","PAM50_genefu","PR.Status","HER2.Final.Status"), cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, cex.dot = 1.5)

plot.new()
mtext("Region by probability (by NormTop)")

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

plot.new()
mtext("Region by Z-score (by NormTop)")

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)


plot.new()
mtext("Region by NormTop (by NormTop)")

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='NormTop', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='NormTop', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)


######################################################
plot.new()
mtext("regionScore by probability")

cisTopicObject <- getRegionsScores(cisTopicObject, method='Probability', scale=TRUE)
#par(mfrow=c(3,5))
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)


cisTopicObject <- GREAT(cisTopicObject, genome='hg19', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
ontologyDotPlot(cisTopicObject, top=5, var.y='name', order.by='Binom_Adjp_BH')


# plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
#      xaxt='n', yaxt='n', xlab='', ylab='')
# topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
#text(1,4, getwd() , pos=4)

#should notify snakemake that pdf is output
cisTopicObject <- runUmap(cisTopicObject, target="cell", n_components=3)
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("PAM50","ER.Status","PAM50_genefu","PR.Status","HER2.Final.Status"), cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, cex.dot = 1.5)

plot.new()
mtext("Region by probability (by probability)")

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

plot.new()
mtext("Region by Z-score (by probability)")

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

plot.new()
mtext("Region by NormTop (by probability)")

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='NormTop', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='NormTop', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)


dev.off()

###Restore warnings####
options(warn = oldw)


######################################THIS IS NOT IN OWN MODULE
topicAssigToPatientOut <- args[6]

#### Write out topic assignments to the patients
write.csv(cisTopicObject@selected.model$document_expects, file=topicAssigToPatientOut)
#write.csv(cisTopicObject@selected.model$document_expects, file="/storage/mathelierarea/processed/petear/analysis/test/topicAssigToPatient.csv")



#### Region scores per topic (normalized umap)
RegScrPrtopicOut <- args[7]
write.csv(cisTopicObject@region.data, file=RegScrPrtopicOut)
#write.csv(cisTopicObject@region.data, file="/storage/mathelierarea/processed/petear/analysis/test/RegScrPrtopic.csv")

#### Unnormalized region assignments
RegAssigUnormalOut <- args[8]
write.csv(cisTopicObject@selected.model$topics, file=RegAssigUnormalOut)
#write.csv(cisTopicObject@selected.model$topics, file="/storage/mathelierarea/processed/petear/analysis/test/RegAssigUnormal.csv")bina



##########################################################################################
library(ComplexHeatmap)
#Create empty dataframe:
dfN <- data.frame(matrix(ncol=2, nrow=1))
colnames(dfN) <- c("chrPos","TopicX")

#merge dataframes
for (attr in attributes(cisTopicObject@binarized.cisTopics)$names)
{
    makeDF <- data.frame(cisTopicObject@binarized.cisTopics[attr])
    rowColDF <- tibble::rownames_to_column(makeDF, "chrPos")
    dfN <- merge(dfN, rowColDF, by="chrPos", all=TRUE)
}
dfN$TopicX <- NULL

binaOut <- args[10]
write.csv(dfN, file=binaOut)
#write.csv(dfN, file="/storage/mathelierarea/processed/petear/analysis/test/bina.csv")

metaOut <- args[11]
write.csv(meta, file=metaOut)
#write.csv(meta, file="/storage/mathelierarea/processed/petear/analysis/test/meta.csv")

# make a matrix with hierarchical clustering as cistopicheatmap :
hmat <- scale(cisTopicObject@selected.model$document_expects, center=TRUE, scale=TRUE)

#make rownames as cistopic gives topic number as rowname (E.g. topic 1 is the first row)
row.names(hmat) <- 1 : nrow(hmat)


rowOrder <- hclust(dist(hmat))$order
colOrder <- hclust(dist(t(hmat)))$order

hmat <- hmat[rowOrder, colOrder]

hmatOut <- args[12]
hmatOut <- write.csv(hmat, file=hmatOut)
#write.csv(hmat, file="/storage/mathelierarea/processed/petear/analysis/test/hmat.csv")



#####get which cluster each patient belong to with help from guidohooiveld on Github #####
 
thamat = t(hmat)                #Transpose hmat because complexHeatmaps k-means clustering (which I have to use for getting cluster assignments) only works on rows.

het <- Heatmap(thamat, km=4)     #Draw a heatmap with km clusters.
#HM <- draw(HM)                  #Show the heatmap
r.dend <- row_dend(het)
rcl.list <- row_order(het)

library(magrittr)           #Probably not needed            

ClusterDF <- lapply(names(rcl.list), function(i){
out <- data.frame(Patient = rownames(thamat[rcl.list[[i]],]),
                  Cluster = paste0("Cluster ", i),
                  stringsAsFactors = FALSE)
return(out)
}) %>%  
do.call(rbind, .)

ClusterDFOut <- args[13]
write.csv(ClusterDF, file=ClusterDFOut)
#write.csv(ClusterDF, file="/storage/mathelierarea/processed/petear/analysis/test/ClusterDF.csv")


####Cutree###   
#Cutree on colorder since they are the Donors

#hmat <- cisTopicObject@selected.model$document_expects

#z-score the hmat:
hmat <- scale(cisTopicObject@selected.model$document_expects, center=TRUE, scale=TRUE)

#make rownames as cistopic gives topic number as rowname (E.g. topic 1 is the first row)
row.names(hmat) <- 1 : nrow(hmat)


rowOrder <- hclust(dist(hmat))$order
colOrder <- cutree(hclust(dist(t(hmat))),k = 1:5)

hmat <- hmat[rowOrder, colOrder]

thamat = t(hmat)                #Transpose hmat because complexHeatmaps k-means clustering (which I have to use for getting cluster assignments) only works on rows.

pdf("hmClus5.pdf", height=15, width=10)

Het <- Heatmap(thamat)     #Draw a heatmap with km clusters.
Het <- draw(Het)                  #Show the heatmap
r.dend <- row_dend(Het)
rcl.list <- row_order(Het)


dev.off()

