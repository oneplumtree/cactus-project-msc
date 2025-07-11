#Snippet of code for gene counts that GAWN outputs to remove the redundant isomers
rm(list=ls())
smallgenes <- read.delim("smallQch1-11_annotation_table.tsv")

annotated <- smallgenes[grep("^-", smallgenes$GeneName,  invert = TRUE), ]

sep1 <- strsplit (annotated[,5],split = "_")
annotated$TranscriptName <- sapply (sep1,"[[",2)

annotatednumbers <- annotated[!duplicated(annotated[5]),]

summary (annotatednumbers) 
