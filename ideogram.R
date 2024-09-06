
#install.packages("RIdeogram")

require(RIdeogram)
source("asynt.R")

### Small Ideogram ----

smallcactus_karyotype <- read.table("QSsmallchromsizes.txt", sep = "\t", header = T, stringsAsFactors = F)

sep1 <- strsplit (smallcactus_karyotype[,1],split = "_")
smallcactus_karyotype$Chr <- sapply (sep1,"[[",3)

ideogram(smallcactus_karyotype, output = "smallchromosome.svg")


gene_density <- GFFex(input = "singlelineWholeQSjour.gff3", karyotype = "smallchromosomesizes.txt", feature = "gene", window = 1000000)
#gene_density <- GFFex(input = "smallgmapJouralkaloidsrenamedchrome.gff3", karyotype = "smallchromosomesizes.txt", feature = "gene", window = 1000000)
alkaloids1 <- read.table("alkaloidideogramwhole.txt", sep = "\t", header = T, stringsAsFactors = F)
alkaloids1 <- alkaloids1[, c(1:5,7)]

alkaloids1$Chr <- as.character(alkaloids1$Chr)

ideogram(karyotype = smallcactus_karyotype, overlaid = gene_density, label = alkaloids1, label_type = "marker", output = "overlaidsmallchromosome.svg")

ideogram(karyotype = smallcactus_karyotype, label = alkaloids1, label_type = "marker", output = "alkaloidsmallchromosome.svg")
convertSVG("alkaloidsmallchromosome.svg", device = "png")

### Big ideogram ----

rawgene_density <- read.table("BiggmapJourmapping.gff3", sep = "\t", header = F, stringsAsFactors = F)
scaffoldrename <- read.table("CondensedScaffolds.txt")

rawgene_density$Chr <- scaffoldrename$V2[match(rawgene_density$V1, scaffoldrename$V1)]
rawgene_density$Chr[is.na(rawgene_density$Chr)] <- "Unscaffolded"

#Fix the coordinates
rawgene_density$adjust <- scaffoldrename$V3[match(rawgene_density$V1, scaffoldrename$V1)]
rawgene_density$Start<- rowSums(rawgene_density[,c("V4", "adjust")])
rawgene_density$End<- rowSums(rawgene_density[,c("V5", "adjust")])

#Remove Unscaffolded and rename
sorted <- rawgene_density[(!(rawgene_density$Chr == "Unscaffolded")),]
split1 <- strsplit (sorted[,9],"Name=")
split2 <- strsplit (sapply (split1,"[[",2),split = ";")
sorted$V10 <- sapply (split2,"[[",1)

sorted2 <- sorted[, c(10,2,3,13,14,6,7,8,15)]
#write.table(sorted2,"Biggene_density.gff3", col.names= F, row.names = F, sep = "\t", quote=F)


toMatch <- c("TRINITY_DN10959_c0_g1_i1",
             "TRINITY_DN10959_c0_g1_i4",
             "TRINITY_DN26_c0_g1_i5",
             "TRINITY_DN9001_c0_g1_i1",
             "TRINITY_DN15955_c0_g2_i1",
             "TRINITY_DN983_c*"
)

matches <- sorted2[grep(paste(toMatch,collapse="|"), 
                        sorted2$V10),]

matches$Type <- "alkaloid" 
matches$Shape <- "circle"
matches$color <- "6a3d9a"

matches2 <- matches[((matches$V3 == "gene")),]

bigalkaloid <- matches2[, c(10,12,1,4,5,11)]
ideogram(karyotype = bigcactus_karyotype, overlaid = biggene_density, label = bigalkaloid, label_type = "marker", output = "bigchromosomemarkers.svg" )
convertSVG("bigchromosomemarkers.svg", device = "png")




bigcactus_karyotype <- read.table("bigScaffoldsizes.txt", sep = "\t", header = T, stringsAsFactors = F)
biggene_density <- GFFex(input = "Biggene_density.gff3", karyotype = "bigScaffoldsizes.txt", feature = "gene", window = 1000000)
ideogram(karyotype = bigcactus_karyotype, overlaid = biggene_density, output = "bigchromosome.svg" )
convertSVG("bigchromosome.svg", device = "png", output = "bigchromosome.png")


### Combined Karyotype ----
cactus_karyotype <- read.table("QSchromsizes.txt", sep = "\t", header = T, stringsAsFactors = F)

sep1 <- strsplit (cactus_karyotype[,1],split = "_")
cactus_karyotype$Chr <- sapply (sep1,"[[",3)

renaming <- read.table("renamekaryotype.txt", header = T, stringsAsFactors = F)

cactus_karyotype$Chr <- renaming$new_karyotype[match(cactus_karyotype$Chr, renaming$old_karyotype)]

library(dplyr)

cactus_karyotype <- cactus_karyotype %>% arrange(Chr)

#ideogram(cactus_karyotype, output = "wholechromosome.svg")
#convertSVG("wholechromosome.svg", device = "png")

#gene_density <- GFFex(input = "singlelineWholeQSjour.gff3", karyotype = "QSchromsizes.txt", feature = "gene", window = 1000000)

gene_densitysmall <- GFFex(input = "singlelinesmallQSjour.gff3", karyotype = "QSchromsizes.txt", feature = "gene", window = 1000000)
gene_densitybig <- GFFex(input = "singlelinebigQSjour.gff3", karyotype = "QSchromsizes.txt", feature = "gene", window = 1000000)

gene_density <- rbind(gene_densitysmall, gene_densitybig)

sep1 <- strsplit (gene_density[,1],split = "_")
gene_density$Chr <- sapply (sep1,"[[",3)

gene_density$Chr <- renaming$new_karyotype[match(gene_density$Chr, renaming$old_karyotype)]


alkaloids1 <- read.table("alkaloidideogramwholenew.txt", sep = "\t", header = T, stringsAsFactors = F)
alkaloids1 <- alkaloids1[, c(1:5,7)]


alkaloids1$Chr <- renaming$new_karyotype[match(alkaloids1$Chr, renaming$old_karyotype)]
alkaloids1$Chr <- as.character(alkaloids1$Chr)


ideogram(karyotype = cactus_karyotype, overlaid = gene_density, label = alkaloids1, label_type = "marker", output = "totalchromosome.svg", Lx = 20, Ly = 25 )
convertSVG("totalchromosome.svg", device = "png")

plot (as.numeric(gene_densitysmall$Value)~as.factor(gene_densitysmall$Chr))
plot (as.numeric(gene_densitybig$Value)~as.factor(gene_densitybig$Chr))
summary(gene_densitysmall)
summary(gene_densitybig)
#write.table(totalcactus_karyotype, file = "totalcactus_karyotype.txt", append = FALSE, quote = FALSE, sep = " ",
 #           eol = "\n", na = "NA", dec = ".", row.names = F,
 #           col.names = TRUE)



##Synteny ----
source("asynt.R")

#install.packages("RIdeogram")
#require("XML")
library(RIdeogram)
# If we have multiple scaffolds making up a chromosome,
# we can string them together, either automatically or manually

#import alignments (note that there are also options import.blast and import.nucmer if you used those tools)
alignments <- import.paf("smallvslarge_wholegenome.paf")

alignments2 <- subset(alignments, Rlen >= 200 & Qlen >= 200)
alignments2 <- subset(alignments2, ref_data$seq_len >= 1000000 & query_data$seq_len >= 1000000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"

subsetsynteny2 <- synteny[grep("_scaffold", synteny$query),]

scaffoldrename <- read.table("CondensedScaffolds.txt")

subsetsynteny2$adjust <- scaffoldrename$V3[match(subsetsynteny2$query, scaffoldrename$V1)]
subsetsynteny2$qStart<- rowSums(subsetsynteny2[,c("Qstart", "adjust")])
subsetsynteny2$qEnd<- rowSums(subsetsynteny2[,c("Qend", "adjust")])

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",2)

sep1 <- strsplit (subsetsynteny2[,1],split = "scaffold")
subsetsynteny2$reference <- sapply (sep1,"[[",2)


subsetsynteny3 <- subsetsynteny2[, c(1:4,9,10,7)]
subsetsynteny3$query <- scaffoldrename$V4[match(subsetsynteny3[,4], scaffoldrename[,1])]


#sep2 <- strsplit (subsetsynteny2[,4],split = "_")
#subsetsynteny2$query <- sapply (sep2,"[[",2)

#sep2 <- strsplit (subsetsynteny2[,4],split = "scaffold")
#subsetsynteny2$query <- sapply (sep2,"[[",2)





#sep1 <- strsplit (subsetsynteny3[,1],split = "_")

#subsetsynteny3$reference <- sapply (sep1,"[[",2)
#sep1 <- strsplit (subsetsynteny3[,1],split = "scaffold")


subsetsynteny3$reference <- as.numeric(subsetsynteny3[,1])
subsetsynteny3$query <- as.numeric(subsetsynteny3[,4])
str(subsetsynteny3)

subsetsynteny3$fill[subsetsynteny3$reference == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11] <- "a71f33"


dualkaryotype22 <- read.table("dualkaryotype1_22.txt", header = T)


ideogram(karyotype = dualkaryotype22, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")

###if there is a temp_table error, it's because the species column is uppercase and not lowercase
dualkaryotype <- read.table("dualkaryotype.txt", header = T)
dualkaryotype50 <- read.table("dualkaryotype1_50.txt", header = T)
table(dualkaryotype50$species)
ideogram(karyotype = dualkaryotype50, synteny = subsetsynteny2)



#Remove Unscaffolded and rename
sorted <- rawgene_density[(!(rawgene_density$Chr == "Unscaffolded")),]
split1 <- strsplit (sorted[,9],"Name=")
split2 <- strsplit (sapply (split1,"[[",2),split = ";")
sorted$V10 <- sapply (split2,"[[",1)


###Small to Beets Synteny ----

alignments <- import.paf("mm2asm20SmallPeyotebeets.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
alignments2 <- subset(alignments2, ref_data$seq_len >= 300000 & query_data$seq_len >= 300000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",2)

sep1 <- strsplit (subsetsynteny2[,1],split = "scaffold")
subsetsynteny2$reference <- sapply (sep1,"[[",2)


dualkaryotype <- read.table("beetsdualkaryotype.txt", header = T)


scaffoldrename <- read.table("BeetsCondensedScaffolds.txt")

subsetsynteny2$query <- scaffoldrename$V2[match(subsetsynteny2[,4], scaffoldrename[,1])]
#sep1 <- strsplit (subsetsynteny3[,1],split = "_")

subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")

subsetsynteny3$fill[subsetsynteny3$reference == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11] <- "a71f33"

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")

###Beets to Small Peyote Synteny ----

alignments <- import.paf("mm2asm20RefBeetssmallPeyote.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
alignments2 <- subset(alignments2, ref_data$seq_len >= 300000 & query_data$seq_len >= 300000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
#subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (synteny[,4],split = "_")
synteny$query <- sapply (sep1,"[[",2)

sep1 <- strsplit (synteny[,4],split = "scaffold")
synteny$query <- sapply (sep1,"[[",2)


dualkaryotype <- read.table("beetsdualkaryotypeinverse.txt", header = T)

scaffoldrename <- read.table("BeetsCondensedScaffolds.txt")

subsetsynteny2 <- synteny

subsetsynteny2$reference <- scaffoldrename$V2[match(synteny[,1], scaffoldrename[,1])]
#sep1 <- strsplit (subsetsynteny3[,1],split = "_")


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, reference > 0)

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")

subsetsynteny3$fill[subsetsynteny3$reference == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11] <- "a71f33"

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")

###Small to Oleraceae Synteny ----

alignments <- import.paf("mm20PeyotesmallQueryOleraceae.paf")

alignments2 <- subset(alignments, Rlen >= 200 & Qlen >= 200)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",2)

sep1 <- strsplit (subsetsynteny2[,1],split = "scaffold")
subsetsynteny2$reference <- sapply (sep1,"[[",2)


dualkaryotype <- read.table("dualkaryotypeoleraceae.txt", header = T)
condensedrename <- read.table("OleraceaeCondensedScaffolds.txt", header = T)

#subsetsynteny2$query <- dualkaryotype$Chr[match(subsetsynteny2[,4], dualkaryotype[,])]
#sep1 <- strsplit (subsetsynteny3[,1],split = "_")

subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill[subsetsynteny3$reference == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11] <- "a71f33"

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")

###Oleraceae to Small Peyote Synteny ----

alignments <- import.paf("mm20RefOleraceaePeyotesmall.paf")

alignments2 <- subset(alignments, Rlen >= 700 & Qlen >= 700)

alignments2 <- subset(alignments2, ref_data$seq_len >= 1000000 & query_data$seq_len >= 1000000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
#subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (synteny[,4],split = "_")
synteny$query <- sapply (sep1,"[[",2)

sep1 <- strsplit (synteny[,4],split = "scaffold")
synteny$query <- sapply (sep1,"[[",2)


dualkaryotype <- read.table("dualkaryotypeoleraceaeinverse.txt", header = T)

scaffoldrename <- read.table("OleraceaeCondensedScaffolds.txt")

subsetsynteny2 <- synteny

subsetsynteny2$reference <- scaffoldrename$V1[match(synteny[,1], scaffoldrename[,2])]
#sep1 <- strsplit (subsetsynteny3[,1],split = "_")


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, reference > 0)


subsetsynteny3$fill[subsetsynteny3$query == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$query == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$query == 11] <- "a71f33"

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")


###Small to Large QS Synteny ----

alignments <- import.paf("mm20RefQSlargesmalljour.paf")

alignments2 <- subset(alignments, Rlen >= 10000 & Qlen >= 10000)
ref_data <- import.genome(fai_file="smallQch1-11.fasta.fai")
query_data <- import.genome(fai_file="LargeQch1-11.fasta.fai")
alignments2 <- subset(alignments2, ref_data$seq_len >= 1000000 & query_data$seq_len >= 1000000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

sep1 <- strsplit (subsetsynteny2[,4],split = "_")
subsetsynteny2$query <- sapply (sep1,"[[",3)


dualkaryotype <- read.table("QSsmallLargedualkaryotypeinverse.txt", header = T)

subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])

renameto <- c(1:11)
currentname <- c(12:22)
condensedrename <- data.frame(renameto,currentname)
subsetsynteny2$reference <- condensedrename$renameto[match(subsetsynteny2[,1], condensedrename[,2])]

str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)


#ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
#convertSVG("chromosome.svg", device = "png")

subsetsynteny3$fill[subsetsynteny3$reference == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11] <- "a71f33"

chrome1 <- subset(subsetsynteny3, reference == 1)

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")


###QS Small to Oleraceae Synteny ----

alignments <- import.paf("mm2RefPeyoteQSsmallOleraceae.paf")

alignments2 <- subset(alignments, Rlen >= 800 & Qlen >= 800)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("dualkaryotypeoleraceaeQS.txt", header = T)
condensedrename <- read.table("OleraceaeCondensedScaffolds.txt", header = T)

subsetsynteny2$query <- condensedrename$X1[match(subsetsynteny2[,4], condensedrename[,2])]
#sep1 <- strsplit (subsetsynteny3[,1],split = "_")

subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill[subsetsynteny3$reference == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11] <- "a71f33"

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")


###QS Small to Dragonfruit Synteny ----

alignments <- import.paf("mm20reorderedQSsmallpitaya.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSsmalldragonfruitkaryotypeshortened.txt", header = T)
condensedrename <- read.table ("renamepitaya.txt")
#subsetsynteny2$query <- gsub('\\D','', subsetsynteny2$query)
#subsetsynteny2$query <-gsub('0','', subsetsynteny2$query)
subsetsynteny2$query <- condensedrename$V2[match(subsetsynteny2[,4], condensedrename[,1])]
subsetsynteny2$reference <- condensedrename$V3[match(subsetsynteny2[,1], condensedrename[,2])]
#sep1 <- strsplit (subsetsynteny3[,1],split = "_")

subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & (subsetsynteny3$query == 9 | subsetsynteny3$query == 8) ] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 9] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]


lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}


sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 8),3])

chromo <- data.frame(subsetsynteny2[,c(1,4)])
table(chromo)
hist(table(chromo))

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")

###QS Dragonfruit to Small Inverse ----

alignments <- import.paf("mm20pitayareorderedQS.paf")

alignments2 <- subset(alignments, Rlen >= 200 & Qlen >= 200)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$query),]

sep1 <- strsplit (subsetsynteny2[,4],split = "_")
subsetsynteny2$query <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSsmalldragonfruitkaryotypeshortenedinverse.txt", header = T)

condensedrename <- read.table ("renamepitaya.txt")
subsetsynteny2$reference <- condensedrename$V2[match(subsetsynteny2[,1], condensedrename[,1])]

df <- read.table("QSrenamebig.txt", header = T)
subsetsynteny2$query <- df$new_small[match(subsetsynteny2[,4], df[,4])]

subsetsynteny2$query <- condensedrename$V4[match(subsetsynteny2[,4], condensedrename[,2])]


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & (subsetsynteny3$query == 9 | subsetsynteny3$query == 8) ] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 9] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]


lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}


sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 8),3])



ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")


chromo <- data.frame(subsetsynteny3[,c(1,4)])
table(chromo)


###QS small and big ----


alignments <- import.paf("mm20genesRefsmallQSbig.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSkaryotypeshortened.txt", header = T)

sep1 <- strsplit (subsetsynteny2[,4],split = "_")
subsetsynteny2$query <- sapply (sep1,"[[",3)

#df <- data.frame(new = 1:11, old = 12:22)
df <- read.table("QSrenamebig.txt", header = T)
subsetsynteny2$query <- df$new_big[match(subsetsynteny2[,4], df[,2])]
subsetsynteny2$reference <- df$new_small[match(subsetsynteny2[,1], df[,4])]


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)
#subsetsynteny3 <- subsetsynteny2[subsetsynteny2$query == 12,]


subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & (subsetsynteny3$query == 9 | subsetsynteny3$query == 2)] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]

chromo <- data.frame(subsetsynteny3[,c(1,4)])

table(chromo)


lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}



sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 2),3])



ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")

###QS big and small  ----


alignments <- import.paf("mm20ReflargeQSsmall.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSkaryotypeshortenedinverse.txt", header = T)

sep1 <- strsplit (subsetsynteny2[,4],split = "_")
subsetsynteny2$query <- sapply (sep1,"[[",3)

#df <- data.frame(new = 1:11, old = 12:22)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)
#subsetsynteny3 <- subsetsynteny2[subsetsynteny2$query == 12,]


subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & subsetsynteny3$query == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]


df <- read.table("QSrenamebasic.txt", header = T)
subsetsynteny2$reference <- df$new_big[match(subsetsynteny2[,1], df[,2])]
subsetsynteny2$query <- df$new_small[match(subsetsynteny2[,4], df[,4])]


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

chromo <- data.frame(subsetsynteny2[,c(1,4)])
table(chromo)


lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}


sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 6),3])

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == (j)),3])
  print(savesum)
}


ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")

#Dragonfruit to Large -----

alignments <- import.paf("mm20refDragonfruitBigQS.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$query),]

sep1 <- strsplit (subsetsynteny2[,4],split = "_")
subsetsynteny2$query <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSdragonfruitlargekaryotype.txt", header = T)

condensedrename <- read.table ("renamepitaya.txt")
subsetsynteny2$reference <- condensedrename$V2[match(subsetsynteny2[,1], condensedrename[,1])]

df <- read.table("QSrenamebig.txt", header = T)
subsetsynteny2$query <- df$df_big[match(subsetsynteny2[,4], df[,2])]

#subsetsynteny2$query <- condensedrename$V4[match(subsetsynteny2[,4], condensedrename[,2])]


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & subsetsynteny3$query == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")


chromo <- data.frame(subsetsynteny3[,c(1,4)])
table(chromo)

lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== 11 & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== 11 & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}


sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 8),3])


#Large Dragonfruit -----

alignments <- import.paf("mm20refQSlargedragonfruit.paf")


alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

condensedrename <- read.table ("renamepitaya.txt")
subsetsynteny2$query <- condensedrename$V2[match(subsetsynteny2[,4], condensedrename[,1])]

df <- read.table("QSrenamebig.txt", header = T)
subsetsynteny2$reference <- df$df_big[match(subsetsynteny2[,1], df[,2])]

subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & subsetsynteny3$query == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]

chromo <- data.frame(subsetsynteny2[,c(1,4)])
table(chromo)

lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)


for (j in 1:11){
print(j)
savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])

print(savesum)
print(oppositesum)
}


n=numeric(0)
list_mat <- list()
for (i in 1:11){
  n[i]=5^i
  m=numeric(0)
  m=matrix(data=0,nrow=n[i],ncol=n[i])
  
  
  for (j in 1:11){
    for (k in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (k)),3])
  print(savesum)
    }
  }
  list_mat[[i]] <- savesum #Holding Matrix
}


sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 8),3])


dualkaryotype <- read.table("QSdragonfruitlargekaryotypeinverse.txt", header = T)


ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")