
rm(list=ls())
smallgenes <- read.delim("smallQch1-11_annotation_table.tsv")

annotated <- smallgenes[grep("^-", smallgenes$GeneName,  invert = TRUE), ]

sep1 <- strsplit (annotated[,5],split = "_")
annotated$TranscriptName <- sapply (sep1,"[[",2)

annotatednumbers <- annotated[!duplicated(annotated[5]),]