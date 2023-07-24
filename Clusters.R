library(ComplexHeatmap)
suppressWarnings(library("cisTopic"))

args <- commandArgs(trailingOnly = TRUE)
argCounter = 1
#Read in bina
#bina <- args[1]



#read in CTO 
cto <- readRDS(args[argCounter])
argCounter = argCounter + 1

hmat <- scale(cto@selected.model$document_expects, center=TRUE, scale=TRUE)
row.names(hmat) <- 1 : nrow(hmat)
rowOrder <- hclust(dist(hmat))$order
colOrder <- hclust(dist(t(hmat)))$order

hmat <- hmat[rowOrder, colOrder]
thamat = t(hmat)                #Transpose hmat because complexHeatmaps k-means clustering (which I have to use for getting cluster assignments) only works on rows.


ClusteredHeatmapsPDF <- args[argCounter]
argCounter = argCounter + 1

pdfNameUMAP <- ClusteredHeatmapsPDF

pdf(pdfNameUMAP, width=10, height=10)
i = 1
while ( i < 10) 
{
i <- i + 1
title_text <- paste("Kmeans clustering = ", i, sep="")
Het <- Heatmap(thamat,show_row_names = FALSE, km=i)     #Draw a heatmap with km clusters.
grid.text(title_text, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))

Het <- draw(Het)                  #Show the heatmap
r.dend <- row_dend(Het)
rcl.list <- row_order(Het)
}
dev.off()


