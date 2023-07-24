######################################################
#####Umap and further cisTopic
######################################################
suppressWarnings(library("cisTopic"))

args <- commandArgs(trailingOnly = TRUE)
#Argcounter is used to keep track of the arguments in case more are added
argCounter = 1

# Load the data from args
ctoIn <- args[argCounter]
argCounter = argCounter + 1


cisTopicObject <- readRDS(ctoIn)


UmapPDF <- args[argCounter]
argCounter = argCounter + 1

pdfNameUmap <- UmapPDF

pdf(pdfNameUmap)

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
