
options(stringsAsFactors = F)

check_overlap_length <- function (to_match,all_start,all_stop,query_start,query_stop){
		
	c1 <- (all_start > query_start & all_start < query_stop) #all_start falls within qstart + qstop
	c2 <- (all_start < query_start & all_start > query_stop) #all_start falls within qstart + qstop [reversed]
	c3 <- (all_stop > query_start & all_stop < query_stop)   #all_stop falls within qstart + qstop
	c4 <- (all_stop < query_start & all_stop > query_stop)   #all_stop falls within qstart + qstop [reversed]
	
	c5 <- (all_start < query_start & all_stop > query_start) #query start falls within all_start and all_stop
	c6 <- (all_start > query_start & all_stop < query_start) #query start falls within all_start and all_stop [reversed]
	c7 <- (all_start < query_stop & all_stop > query_stop)   #query stop falls within all_start and all_stop
	c8 <- (all_start > query_stop & all_stop < query_stop)   #query stop falls within all_start and all_stop [reversed]
		
	test1 <- to_match[(c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8),]

	sub_astart <- all_start[(c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8)]
	sub_astop <- all_stop[(c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8)]

	query_range <- query_start:query_stop
	

	if (nrow (test1) > 0){
		test1$overlap <- NA

		for (yyy in 1:nrow(test1)){
			vec1 <- sub_astart[yyy]:sub_astop[yyy]
			
			test1$overlap[yyy] <- sum (vec1 %in% query_range)
					
		}
		
	}
	
	test1
}


##Original Code ------
jour1 <- read.table ("Jour_wholejour.gtf",sep = "\t")
will1 <- read.table ("Williamsii_wholejour.gtf",sep = "\t")

jour1b <- jour1[grep ("TPM", jour1[,9]),]
will1b <- will1[grep ("TPM", will1[,9]),]

jour1b <- jour1b[grep ("PGA", jour1b[,1]),]
will1b <- will1b[grep ("PGA", will1b[,1]),]



sep1 <- strsplit (jour1b[,9],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",5),split = "TPM ")
jour1b$tpm <- sapply (sep2,"[[",2)


sep1 <- strsplit (will1b[,9],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",5),split = "TPM ")
will1b$tpm <- sapply (sep2,"[[",2)


chromsplit1 <- strsplit (will1b[,1],split = "__")
chromsplit2 <- strsplit (sapply (chromsplit1,"[[",1),"scaffold")
will1b$chromnumeric <- as.numeric (sapply(chromsplit2,"[[",2))

chromsplit1 <- strsplit (jour1b[,1],split = "__")
chromsplit2 <- strsplit (sapply (chromsplit1,"[[",1),"scaffold")
jour1b$chromnumeric <- as.numeric (sapply(chromsplit2,"[[",2))


#reorder
will1b <- will1b[order(will1b$chromnumeric,will1b[,4]),]
jour1b <- jour1b[order(jour1b$chromnumeric,jour1b[,4]),]


willmean1 <- mean (as.numeric (will1b$tpm))
jourmean1 <- mean (as.numeric (jour1b$tpm))

par (mfcol = c (2,2), mar = c(1,1,1,1))
plot (will1b[,10], ylim = c (0,250))
arrows(-1000,willmean1,1000000,willmean1,col = "red")
plot (jour1b[,10], ylim = c (0,250))
arrows(-1000,jourmean1,1000000,jourmean1,col = "red")


dim (will1b)
dim (jour1b)


will1b$ID <- paste (will1b$chromnumeric,will1b[,4],sep = "__")

jour1b$ID <- paste (jour1b$chromnumeric,jour1b[,4],sep = "__")

merg1 <- merge (will1b,jour1b,by.x = "ID",by.y = "ID")


par (mfcol = c(1,2))
plot (merg1$tpm.x,merg1$tpm.y,xlim = c (0,300),ylim = c (0,300))
arrows(-1000,-1000,1000000,1000000,col = "red")

#sub_jour1b <- jour1b[as.numeric (jour1b$tpm) > 1000,]

jour1b$will1b_tpm <- NA

for (i in 1:nrow(jour1b)){
	
	sub_test <- will1b[will1b[,1] == jour1b[i,1],]	
	
	test1 <- check_overlap_length (sub_test,as.numeric (sub_test[,4]),as.numeric (sub_test[,5]),as.numeric (jour1b[i,4]),as.numeric (jour1b[i,5]))
	
	jour1b$will1b_tpm[i] <- mean (as.numeric (test1$tpm))
	
	
}




par (mfcol = c(1,1), mar = c(5,5,3,3))
plot (log(as.numeric(smalljour$tpm)),log(as.numeric(smalljour$will1b_tpm)),xlim = c(-1,6),ylim = c(-1,6),xlab = "Jourdaniana Small Genome Log TPM ",ylab = "Williamsii Small Genome Log TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")

plot (log(as.numeric(bigjour$tpm)),log(as.numeric(bigjour$will1b_tpm)),xlim = c(-1,6),ylim = c(-1,6),xlab = "Jourdaniana Big Genome Log TPM ",ylab = "Williamsii Big Genome Log TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")

plot (smalljour$tpm,smalljour$will1b_tpm,xlim = c(0,100),ylim = c(0,100),xlab = "Jourdaniana Small Genome TPM ",ylab = "Williamsii Small Genome TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")

plot (bigjour$tpm,bigjour$will1b_tpm,xlim = c(0,200),ylim = c(0,200),xlab = "Jourdaniana Big Genome TPM ",ylab = "Williamsii Big Genome TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")

###subset with dots 
plot (smalljour$tpm,smalljour$will1b_tpm,xlim = c(0,100),ylim = c(0,100),xlab = "Jourdaniana Small Genome TPM ",ylab = "Williamsii Small Genome TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")

plot (bigjour$tpm,bigjour$will1b_tpm,xlim = c(0,200),ylim = c(0,200),xlab = "Jourdaniana TPM Reads Mapped to Large Compartment",ylab = "Williamsii TPM Reads Mapped to Large Compartment")
arrows(-1000,-1000,1000000,1000000,col = "red")


library(dplyr)

matching <- read.table("startcoordinates.txt", header = F)

smallmatch <- smalljour[smalljour$V4 %in% matching$V3,]
bigmatch <- bigjour[bigjour$V4 %in% matching$V3,]

smallmatch$enzyme <- matching$V1[match(smallmatch$V4, matching$V3)]
bigmatch$enzyme <- matching$V1[match(bigmatch$V4, matching$V3)]

plot (log(as.numeric(smalljour$tpm)),log(as.numeric(smalljour$will1b_tpm)),xlim = c(0,6),ylim = c(0,6),xlab = "Jourdaniana Small Genome Log TPM ",ylab = "Williamsii Small Genome Log TPM", col=ifelse(smalljour$tpm %in% smallmatch$tpm, 'red', 'grey'))
points(log(as.numeric(smallmatch$tpm)), log(as.numeric(smallmatch$will1b_tpm)), col="red")
arrows(-1000,-1000,1000000,1000000,col = "red")

###This is the coloured points code
col1 <- as.factor(bigmatch$enzyme)
plot(log(as.numeric(bigjour$tpm)),log(as.numeric(bigjour$will1b_tpm)),xlim = c(0,6),ylim = c(0,6),xlab = "Jourdaniana Big Genome Log TPM ",ylab = "Williamsii Big Genome Log TPM", col=ifelse(bigjour$tpm %in% bigmatch$tpm, 'red', 'grey'))
points(log(as.numeric(bigmatch$tpm)), log(as.numeric(bigmatch$will1b_tpm)), col=col1)
arrows(-1000,-1000,1000000,1000000,col = "red")
#legend("topright", inset=c(-0.2,0), legend= c("LwCYP76AD94", "LwOMT2", "LwOMT6", "LwTyDC1"), col = 1:4, cex = 0.8, pch = 1)

plot (smalljour$tpm,smalljour$will1b_tpm,xlim = c(0,100),ylim = c(0,100),xlab = "Jourdaniana Small Genome TPM ",ylab = "Williamsii Small Genome TPM", col=ifelse(smalljour$tpm %in% smallmatch$tpm, 'red', 'grey'))
arrows(-1000,-1000,1000000,1000000,col = "red")
points(smallmatch$tpm, smallmatch$will1b_tpm, col=col1)

plot (bigjour$tpm,bigjour$will1b_tpm,xlim = c(0,200),ylim = c(0,200),xlab = "Jourdaniana Big Genome TPM ",ylab = "Williamsii Big Genome TPM", col=ifelse(bigjour$tpm %in% bigmatch$tpm, 'red', 'grey'))
arrows(-1000,-1000,1000000,1000000,col = "red")
points((bigmatch$tpm), (bigmatch$will1b_tpm), col=col1)



### QS Data =========


jour1 <- read.table ("readsQSjourwhole.gtf",sep = "\t")
will1 <- read.table ("JourWholeQSWilliamsiireads.gtf",sep = "\t")

jour1b <- jour1[grep ("TPM", jour1[,9]),]
will1b <- will1[grep ("TPM", will1[,9]),]

#jour1b <- jour1b[grep ("PGA", jour1b[,1]),]
#will1b <- will1b[grep ("PGA", will1b[,1]),]



sep1 <- strsplit (jour1b[,9],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",5),split = "TPM ")
jour1b$tpm <- sapply (sep2,"[[",2)


sep1 <- strsplit (will1b[,9],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",5),split = "TPM ")
will1b$tpm <- sapply (sep2,"[[",2)


chromsplit1 <- strsplit (will1b[,1],split = "_")
chromsplit2 <- strsplit (sapply (chromsplit1,"[[",3),"scaffold")
will1b$chromnumeric <- as.numeric (sapply(chromsplit2,"[[",1))

chromsplit1 <- strsplit (jour1b[,1],split = "_")
chromsplit2 <- strsplit (sapply (chromsplit1,"[[",3),"scaffold")
jour1b$chromnumeric <- as.numeric (sapply(chromsplit2,"[[",1))


#reorder
will1b <- will1b[order(will1b$chromnumeric,will1b[,4]),]
jour1b <- jour1b[order(jour1b$chromnumeric,jour1b[,4]),]


willmean1 <- mean (as.numeric (will1b$tpm))
jourmean1 <- mean (as.numeric (jour1b$tpm))

plot (will1b[,10], ylim = c (0,250))
arrows(-1000,willmean1,1000000,willmean1,col = "red")
plot (jour1b[,10], ylim = c (0,250))
arrows(-1000,jourmean1,1000000,jourmean1,col = "red")


dim (will1b)
dim (jour1b)


will1b$ID <- paste (will1b$chromnumeric,will1b[,4],sep = "__")

jour1b$ID <- paste (jour1b$chromnumeric,jour1b[,4],sep = "__")

merg1 <- merge (will1b,jour1b,by.x = "ID",by.y = "ID")


plot (merg1$tpm.x,merg1$tpm.y,xlim = c (0,300),ylim = c (0,300))
arrows(-1000,-1000,1000000,1000000,col = "red")

#sub_jour1b <- jour1b[as.numeric (jour1b$tpm) > 1000,]

jour1b$will1b_tpm <- NA

for (i in 1:nrow(jour1b)){
  
  sub_test <- will1b[will1b[,1] == jour1b[i,1],]	
  
  test1 <- check_overlap_length (sub_test,as.numeric (sub_test[,4]),as.numeric (sub_test[,5]),as.numeric (jour1b[i,4]),as.numeric (jour1b[i,5]))
  
  jour1b$will1b_tpm[i] <- mean (as.numeric (test1$tpm))
  
  
}


smalljour <- jour1b[jour1b$chromnumeric < 12 ,]
bigjour <- jour1b[jour1b$chromnumeric > 13 ,]

smalljour$will1b_tpm[is.na(smalljour$will1b_tpm)] <- "0"
bigjour$will1b_tpm[is.na(bigjour$will1b_tpm)] <- "0"
bigjour$tpm[is.na(bigjour$tpm)] <- "0"


library(dplyr)

matching <- read.table("QSstartcoordinates.txt", header = F)

smallmatch <- smalljour[smalljour$V4 %in% matching$V3,]
bigmatch <- bigjour[bigjour$V4 %in% matching$V3,]

smallmatch$enzyme <- matching$V1[match(smallmatch$V4, matching$V3)]
bigmatch$enzyme <- matching$V2[match(bigmatch$V4, matching$V3)]

###This is the coloured points code

plot (log(as.numeric(smalljour$tpm)),log(as.numeric(smalljour$will1b_tpm)),xlim = c(0,6),ylim = c(0,6),xlab = "Jourdaniana Small Genome Log TPM ",ylab = "Williamsii Small Genome Log TPM", col=ifelse(smalljour$tpm %in% smallmatch$tpm, 'red', 'grey'))
points(log(as.numeric(smallmatch$tpm)), log(as.numeric(smallmatch$will1b_tpm)), col="red")
arrows(-1000,-1000,1000000,1000000,col = "red")

par (mfcol = c (1,1), mar = c(5,5,3,3))
col1 <- as.factor(bigmatch$enzyme)
plot(log(as.numeric(bigjour$tpm)),log(as.numeric(bigjour$will1b_tpm)),xlim = c(-1,6),ylim = c(-1,6),xlab = "Jourdaniana Log TPM: Reads Mapped to Large Compartment",ylab = "Williamsii Log TPM: Reads Mapped to Large Compartment", col=ifelse(bigjour$tpm %in% bigmatch$tpm, 'red', 'grey'))
points(log(as.numeric(bigmatch$tpm)), log(as.numeric(bigmatch$will1b_tpm)), col=col1)
arrows(-80,-81,80,81,col = "red")
arrows(log(0.4),-100,log(0.8),100,col = "orange")
arrows(-10,-0.8,10,0.8,col = "orange")
legend("topright", legend= c("LwCYP76AD94", "LwOMT2", "LwOMT6", "LwTyDC1"), col = 1:4, cex = 0.8, pch = 1)

log(as.numeric(bigjour$tpm))

sd(log(as.numeric(bigmatch$tpm)))
mean (log(as.numeric(bigjour$tpm)))

std <- sd(log(as.numeric(bigmatch$will1b_tpm)))
ave <- mean (log(as.numeric(bigmatch$will1b_tpm)))
      
#'real std'      
plot (bigjour$tpm,bigjour$will1b_tpm,xlim = c(0,150),ylim = c(0,150),xlab = "Jourdaniana TPM: Reads Mapped to Large Compartment",ylab = "Williamsii TPM: Reads Mapped to Large Compartment", col=ifelse(bigjour$tpm %in% bigmatch$tpm, 'red', 'grey'))
arrows(-1000,-1000,1000000,1000000,col = "red")
arrows(-2.225541,-22026.47,2.225541,22026.47,col = "orange")
arrows(-22026.47,-2.225541,22026.47,2.225541,col = "orange")
points((bigmatch$tpm), (bigmatch$will1b_tpm), col=col1)

#fake std
plot (bigjour$tpm,bigjour$will1b_tpm,xlim = c(0,150),ylim = c(0,150),xlab = "Jourdaniana TPM: Reads Mapped to Large Compartment",ylab = "Williamsii TPM: Reads Mapped to Large Compartment", col=ifelse(bigjour$tpm %in% bigmatch$tpm, 'red', 'grey'))
arrows(-1000,-1000,1000000,1000000,col = "red")
arrows(-50,-200,50,200,col = "orange")
arrows(-200,-50,200,50,col = "orange")
points((bigmatch$tpm), (bigmatch$will1b_tpm), col=col1)


sdj <- sd(as.numeric(bigmatch$tpm))
avej<- mean (as.numeric(bigjour$tpm))

std <- sd(as.numeric(bigmatch$will1b_tpm))
ave <- mean (as.numeric(bigmatch$will1b_tpm))


bigjour$tpm <- as.numeric (bigjour$tpm)
bigjour$will1b_tpm <- as.numeric (bigjour$will1b_tpm)

willspecial <- subset(bigjour, tpm < 50 & will1b_tpm > 100)
jourspecial <- subset(bigjour, tpm > 40 & will1b_tpm < 25)

library(ape)
transcriptsdata <- read.gff("singlelineWholeQSjour.gff3", na.strings = c(".", "?"), GFF3 = TRUE)
transcriptsdata[,1] <- as.character(transcriptsdata[,1])

chromsplit1 <- strsplit (transcriptsdata[,1],split = "_")
transcriptsdata$chromnumeric <- as.numeric (sapply(chromsplit1,"[[",3))

transcriptsdatabig <- transcriptsdata[transcriptsdata$chromnumeric > 11,]

transcriptsdatabig <- transcriptsdata[transcriptsdata$type == "gene" ,]



datalist <- NA
numbers <- NA
options(stringsAsFactors = F)
for (i in 1:nrow(jourspecial)){
  
  sub_test <- transcriptsdatabig[transcriptsdatabig[,10] == jourspecial[i,11],]	
  test1 <- check_overlap_length (sub_test,as.numeric (sub_test[,4]),as.numeric (sub_test[,5]),as.numeric (jourspecial[i,4]),as.numeric (jourspecial[i,5]))
  datalist <- rbind(datalist,test1)
}


datalist 

datalist2 <- datalist[-1,]

sep1 <- strsplit (datalist2[,9],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",2),split = "ID=")
sep3 <- strsplit (sapply (sep2,"[[",1),split = "Name=")
datalist2$attributes <- sapply (sep3,"[[",2)




wdatalist <- NA
options(stringsAsFactors = F)
for (i in 1:nrow(willspecial)){
  
  sub_test <- transcriptsdatabig[transcriptsdatabig[,10] == willspecial[i,11],]	
  test1 <- check_overlap_length (sub_test,as.numeric (sub_test[,4]),as.numeric (sub_test[,5]),as.numeric (willspecial[i,4]),as.numeric (willspecial[i,5]))
  
  wdatalist <- rbind(wdatalist,test1)
}
wdatalist 

wdatalist2 <- wdatalist[-1,]

sep1 <- strsplit (wdatalist2[,9],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",2),split = "Name=")
wdatalist2$attributes <- sapply (sep2,"[[",2)
View(wdatalist2$attributes)

trintemp <- read.table ("tempTrinity.txt")
left <- datalist2[!datalist2$attributes %in% trintemp$V1,] 
left$attributes


toget <- read.table("toget.txt")
chrom <- wdatalist2[wdatalist2$attributes %in% toget$V1,] 
jourspecial$chromnumeric <- as.factor(jourspecial$chromnumeric)
aggregate (jourspecial$tpm ~ jourspecial$chromnumeric, FUN = mean)
aggregate (jourspecial$will1b_tpm ~ jourspecial$chromnumeric, FUN = mean)


aggregate (willspecial$tpm ~ willspecial$chromnumeric, FUN = mean)
aggregate (willspecial$will1b_tpm ~ willspecial$chromnumeric, FUN = mean)

# Junk Code ---------------------------------


plot (smalljour$tpm,smalljour$will1b_tpm,xlim = c(0,100),ylim = c(0,100),xlab = "Jourdaniana Small Genome TPM ",ylab = "Williamsii Small Genome TPM", col=ifelse(smalljour$tpm %in% smallmatch$tpm, 'red', 'grey'))
arrows(-1000,-1000,1000000,1000000,col = "red")
points(smallmatch$tpm, smallmatch$will1b_tpm, col=col1)


par (mfcol = c(1,1), mar = c(5,5,3,3))
plot (log(as.numeric(smalljour$tpm)),log(as.numeric(smalljour$will1b_tpm)),xlim = c(-1,6),ylim = c(-1,6),xlab = "Jourdaniana Small Genome Log TPM ",ylab = "Williamsii Small Genome Log TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")

plot (log(as.numeric(bigjour$tpm)),log(as.numeric(bigjour$will1b_tpm)),xlim = c(-1,6),ylim = c(-1,6),xlab = "Jourdaniana Big Genome Log TPM ",ylab = "Williamsii Big Genome Log TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")

plot (smalljour$tpm,smalljour$will1b_tpm,xlim = c(0,100),ylim = c(0,100),xlab = "Jourdaniana Small Genome TPM ",ylab = "Williamsii Small Genome TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")

plot (bigjour$tpm,bigjour$will1b_tpm,xlim = c(0,200),ylim = c(0,200),xlab = "Jourdaniana TPM Reads Mapped to Large Compartment",ylab = "Williamsii TPM Reads Mapped to Large Compartment")
arrows(-1000,-1000,1000000,1000000,col = "red")

###subset with dots 
plot (smalljour$tpm,smalljour$will1b_tpm,xlim = c(0,100),ylim = c(0,100),xlab = "Jourdaniana Small Genome TPM ",ylab = "Williamsii Small Genome TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")

plot (bigjour$tpm,bigjour$will1b_tpm,xlim = c(0,200),ylim = c(0,200),xlab = "Jourdaniana Big Genome TPM ",ylab = "Williamsii Big Genome TPM")
arrows(-1000,-1000,1000000,1000000,col = "red")



smalljour <- jour1b[jour1b$chromnumeric < 12 ,]
bigjour <- jour1b[jour1b$chromnumeric > 13 ,]

smalljour$will1b_tpm[is.na(smalljour$will1b_tpm)] <- "0"
bigjour$will1b_tpm[is.na(bigjour$will1b_tpm)] <- "0"
bigjour$tpm[is.na(bigjour$tpm)] <- "0"


hist((as.numeric(smalljour$will1b_tpm)), 
     main="Counts of Williamsii TPM Small",
     xlab="TPM",)

hist((as.numeric(smalljour$tpm)), 
     main="Counts of Jourdaniana TPM Small",
     xlab="TPM",)

hist((as.numeric(bigjour$will1b_tpm)), 
     main="Counts of Williamsii TPM Large",
     xlab="TPM",)

hist((as.numeric(bigjour$tpm)), 
     main="Counts of Jourdaniana TPM Large",
     xlab="TPM",)

d <- density((as.numeric(bigjour$will1b_tpm)))
plot(d)

d <- density((as.numeric(bigjour$tpm)))
plot(d)