
#Load libraries
library(ComplexHeatmap)
library(dplyr)

#Load data
df <- read.csv(snakemake@input[["py_DF"]], header=TRUE, sep="\t")
meta <- read.table(snakemake@input[["py_meta"]], header=TRUE, sep=",")

# df <- read.csv("py_df.csv", header=TRUE, sep="\t")
# meta <- read.csv("py_df.csv", header=TRUE, sep=",")

meta <- subset(meta, subset = -X)

#The columns are not sorted for df, so sorting both to make sure 
df <- df %>% select(sort(names(df)))
meta <- meta %>% arrange(Donors)

#change NA with string NA
meta <- meta %>% replace(is.na(.), "NA")

