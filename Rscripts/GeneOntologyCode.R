rm(list=ls(all=TRUE))

#load libraries
library(ape)
library(dplyr)
library(tidyr)
library(readr)
#install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(IRanges)
library(ggplot2)

#Isolate the gmap file to only small compartment hits
wholegmap <- read.gff("jourQSwholegmapuncondensed.gff3")

smallonly <- wholegmap[wholegmap$seqid %in% c("HiC_scaffold_1", "HiC_scaffold_2", "HiC_scaffold_3", "HiC_scaffold_4", "HiC_scaffold_5", "HiC_scaffold_6", "HiC_scaffold_7", "HiC_scaffold_8", "HiC_scaffold_9", "HiC_scaffold_10", "HiC_scaffold_11"),]
smallonlygenes <- smallonly[smallonly$type == "gene",]

split1 <- strsplit(smallonlygenes$attributes, ";")
split2 <- strsplit(sapply (split1,"[[",2),split = "Name=")
split3 <- strsplit(sapply (split2,"[[",2),split = "_i")
smallonlygenes$ID <- sapply(split3,"[[",1)

smallonlygenesunique <- smallonlygenes[!duplicated(smallonlygenes$ID),]
#####Finding Expression ####

smallonlygenes$tpm <- 0

##clean up file names
stringtieall <- read.table ("readsQSjourwhole.gtf",sep = "\t")
stringtieall <- stringtieall[stringtieall$V3=="transcript",]

sep1 <- strsplit (stringtieall[,9],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",5), split = "TPM ")
stringtieall$tpm <- sapply (sep2,"[[",2)

sep1 <- strsplit (stringtieall[,9],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",2), split = "transcript_id ")
stringtieall$ID <- sapply (sep2,"[[",2)

#find overlapping ranges
smallgmaprange <- GRanges(seqnames = smallonlygenes$seqid,
                          ranges = IRanges(start = smallonlygenes$start, end = smallonlygenes$end),
                          TPM = smallonlygenes$tpm,
                          ID = smallonlygenes$ID)

stringtierange <- GRanges(seqnames = stringtieall$V1,
                          ranges = IRanges(start = stringtieall$V4, end = stringtieall$V5),
                          TPM = stringtieall$tpm,
                          ID = stringtieall$ID)

overlapping_data <- mergeByOverlaps(smallgmaprange, stringtierange, ignore.strand = TRUE)

# View the resulting data frame with TPM values
smallgenevalues <- as.data.frame(overlapping_data)
smallgenevalues[,17] <- as.numeric(smallgenevalues[,17])

#clean up file to reduce extraneous columns
reducedgenessmall <- smallgenevalues[,c(1:3, 9, 16:17)]

#get an aggregated list of unique hits
smallsumTPMname <- aggregate(smallgenevalues$TPM.1~smallgenevalues[,7],smallgenevalues,sum) ##will collapse the ones on different chr into each other 
smallsumTPM <- aggregate(smallgenevalues$TPM.1~smallgenevalues[,3],smallgenevalues,sum) ###proxy bc the end coor unlikley to be the same
#smallaveTPM <- aggregate(smallgenevalues$TPM.1~smallgenevalues[,7],smallgenevalues,mean)

#hist(smallaveTPM[,2])
#hist(smallsumTPM[,2])


#####Finding GO Terms Overlap#####

gawnGO <- read_tsv("smallQch1-11_annotation_table.tsv") #GO terms generated with GAWN, uniprot data.

#library(GenomicRanges)
#library(IRanges)


smallgmaprange <- GRanges(seqnames = smallonlygenes$seqid,
                          ranges = IRanges(start = smallonlygenes$start, end = smallonlygenes$end),
                          ID = smallonlygenes$ID)

gawnrange <- GRanges(seqnames = gawnGO$ScaffoldName,
                          ranges = IRanges(start = gawnGO$FromPosition, end = gawnGO$ToPosition),
                          ID = gawnGO$GeneName,
                          GO = gawnGO$BiologicalProcess,
                          Accession = gawnGO$GeneGo)

# Find overlaps
overlapping_data <- mergeByOverlaps(smallgmaprange, gawnrange, ignore.strand = TRUE)

# View the resulting data frame with GO Annotations
smallGO <- as.data.frame(overlapping_data)

GOfiltered <- smallGO[smallGO$Accession != "-",]

GOunique <-GOfiltered[!duplicated(GOfiltered$smallgmaprange.ID),]


GOuniqueaccessions <- separate_rows(GOunique, Accession, sep = ";")

GOuniqueaccessions2 <- GOuniqueaccessions[GOuniqueaccessions$Accession != "",]
GOuniqueaccessions2 <- GOuniqueaccessions2[,c(1:4, 7, 16:18)]

GOresults <- GOuniqueaccessions2 %>% group_by(smallgmaprange.seqnames, Accession) %>%
  summarise(count = n(), .groups = "drop")

GOresults2 <- GOresults[GOresults$count > 20,] ##change number if you want to sort
#write.csv(GOresults2, "GOresults.csv")
#write.csv(GOunique, "GoTermsSmallcompartment.csv") #wrote the files to copy and paste them into GO website to look up their annotations https://www.ebi.ac.uk/QuickGO/#

uniprotGO <- read_tsv("GObasket.tsv") #exported this list from the site 

GOresults2$Process <- uniprotGO$Annotation[match(GOresults2$Accession, uniprotGO$Gene)]
GOresults2$Category <- uniprotGO$Category[match(GOresults2$Accession, uniprotGO$Gene)]


GOresults2$smallgmaprange.seqnames <- droplevels(GOresults2$smallgmaprange.seqnames)
GOresults2$smallgmaprange.seqnames <- as.character(GOresults2$smallgmaprange.seqnames)

str(GOresults2)

sep1 <- strsplit (GOresults2$smallgmaprange.seqnames,split = "_")
GOresults2$chr <- sapply (sep1,"[[",3)
GOresults2$chr <- as.numeric(GOresults2$chr)

boxplot(GOresults2$count~GOresults2$chr, xlab = "Chromosome", ylab = "Count")
GOresults2$chr <- as.factor(GOresults2$chr)

cellsonly <- GOresults2[GOresults2$Category == "cellular_component",]

library(ggplot2)
Cellsplot <- ggplot(cellsonly, aes(x=chr, y=count, color = Process)) + 
  geom_point()+
  theme_classic()+
  scale_x_discrete(name="Chromosome")+
  scale_y_continuous(name="Count")


molecular <- GOresults2[GOresults2$Category == "molecular_function",]

ggplot(molecular, aes(x=chr, y=count, color = Process)) + 
  geom_point()+
  theme_classic()+
  scale_x_discrete(name="Chromosome")+
  scale_y_continuous(name="Count")


biological <- GOresults2[GOresults2$Category == "biological_process",]

biological2 <- biological[biological$count > 50,]

ggplot(biological2, aes(x=chr, y=count, color = Process)) + 
  geom_point()+
  theme_classic()+
  scale_x_discrete(name="Chromosome")+
  scale_y_continuous(name="Count")


##### Combine GO with Expression #####


stringtierange <- GRanges(seqnames = stringtieall$V1,
                          ranges = IRanges(start = stringtieall$V4, end = stringtieall$V5),
                          TPM = stringtieall$tpm,
                          ID = stringtieall$ID)

gawnrange <- GRanges(seqnames = gawnGO$ScaffoldName,
                     ranges = IRanges(start = gawnGO$FromPosition, end = gawnGO$ToPosition),
                     ID = gawnGO$GeneName,
                     GO = gawnGO$BiologicalProcess,
                     Accession = gawnGO$GeneGo)

# Find overlaps
overlapping_data <- mergeByOverlaps(stringtierange, gawnrange, ignore.strand = TRUE)

# View the resulting data frame with GO Annotations
stringtieGO <- as.data.frame(overlapping_data)

stringtieGOfiltered <- stringtieGO[stringtieGO$Accession != "-",]
stringtieGOunique <- stringtieGOfiltered[!duplicated(stringtieGOfiltered$stringtierange.ID),]

stringtieGOuniqueaccessions <- separate_rows(stringtieGOunique, Accession, sep = ";")

#stringtieannotations <- separate_rows(stringtieGOunique, GO, sep = ";") ###not doing GO P: bc it doesn't annotate everything


stringtieGOuniqueaccessions2 <- stringtieGOuniqueaccessions[stringtieGOuniqueaccessions$Accession != "",]
stringtieGOuniqueaccessions2 <- stringtieGOuniqueaccessions2[,c(1:4, 6:8, 11,12, 16:20)]

strGOresults <- stringtieGOuniqueaccessions2 %>% group_by(stringtierange.seqnames, Accession, TPM) %>%
  summarise(count = n(), .groups = "drop")

#strGOresults <- strGOresults[strGOresults$count > 20,] ##change number if you want to sort
#write.csv(strGOresults, "strGOresults.csv")
write.csv(stringtieGOunique, "strGoTermsSmallcompartment.csv")

struniprotGO <- read_tsv("strbasket.tsv") ###Need to add headers to the tsv file manually!

strGOresults$Process <- struniprotGO$Annotation[match(strGOresults$Accession, struniprotGO$Gene)]
strGOresults$Category <- struniprotGO$Category[match(strGOresults$Accession, struniprotGO$Gene)]


strGOresults$stringtierange.seqnames <- droplevels(strGOresults$stringtierange.seqnames)
strGOresults$stringtierange.seqnames <- as.character(strGOresults$stringtierange.seqnames)

strGOresults$TPM <- as.numeric(strGOresults$TPM)

str(strGOresults)

#visualize results in plots
sep1 <- strsplit (strGOresults$stringtierange.seqnames,split = "_")
strGOresults$chr <- sapply (sep1,"[[",3)
strGOresults$chr <- as.numeric(strGOresults$chr)

boxplot(strGOresults$TPM~strGOresults$chr, xlab = "Chromosome", ylab = "TPM", main = "Stringtie GO Counts")

strfilteredcount <- strGOresults[strGOresults$TPM > 30,]

strcellsonly <- strfilteredcount[strfilteredcount$Category == "cellular_component" & !is.na(strfilteredcount$Process),]

library(ggplot2)
strGOresults$chr <- as.factor(strGOresults$chr)
ggplot(strcellsonly, aes(x=chr, y=TPM, color = Process)) + 
  geom_point()+
  theme_classic()+
  scale_x_discrete(name="Chromosome")+
  scale_y_continuous(name="TPM")


strmolecular <- strGOresults[strGOresults$Category == "molecular_function" & !is.na(strGOresults$Process) & strGOresults$TPM > 65,]

ggplot(strmolecular, aes(x=chr, y=TPM, color = Process)) + 
  geom_point()+
  theme_classic()+
  scale_x_discrete(name="Chromosome")+
  scale_y_continuous(name="TPM")+
  ggtitle("Chromosome 4 Top Molecular Functions")


strbiological <- strGOresults[strGOresults$Category == "biological_process" & !is.na(strGOresults$Process) & strGOresults$TPM > 50 & strGOresults$chr == "4",]

ggplot(strbiological, aes(x=chr, y=TPM, color = Process)) + 
  geom_point()+
  theme_classic()+
  scale_x_discrete(name="Chromosome")+
  scale_y_continuous(name="TPM")
