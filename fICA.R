library(fastICA)
library("ggplot2")
library("RColorBrewer")
suppressWarnings(library("cisTopic"))

#Argcounter is used to keep track of the arguments in case more are added

ctoIn <- snakemake@input[["ctoIn"]]

cisTopicObject <- readRDS(ctoIn)

fIcaPDF <- snakemake@output[["fICApdf"]]

pdf(fIcaPDF)

#Get the document_expects matrix from the cisTopicObject and transpose it
t_cisTo <- t(cisTopicObject@selected.model$document_expects)

#Run ICA on the transposed matrix
donorMat  <- matrix(rownames(t_cisTo))
fica <- fastICA(t_cisTo, n.comp= 2, row.norm=FALSE, method = "C")

#Use only the S component of the ICA
ficaS <- fica$S

#Add V1 and V2 as column names and donor names as row names
colnames(ficaS) <- c("V1", "V2")
rownames(ficaS) <- c(donorMat)


metaficaS <- meta[rownames(meta) %in% rownames(ficaS), ]

#Merge the meta data with the ICA data
metaficaS <- as.data.frame(metaficaS)
ficaS <- as.data.frame(ficaS)
ficaSMerge <- merge(metaficaS, ficaS, by = 0)

library(RColorBrewer)
library(ggplot2)

################# GGPLOT2 for PAM50 ################
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


dev.off() 
