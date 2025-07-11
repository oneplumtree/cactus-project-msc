#Ginny Li - XCMS Runs - H1-C1 March 2025 Runs

rm(list=ls())

library(xcms)
library(RColorBrewer)
library(readMzXmlData)
library(MsExperiment)

#Use Proteowizard MZConvert to change raw to mzxml
#####Hseries ####
xcmsXML <- list.files("C:/Users/ginlo/OneDrive/Desktop/March18-2025/XMLfiles/Hseries", full.names = T, recursive = T) [c(1:18)]

pd <- data.frame(sample_name = sub(basename(xcmsXML), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE),
                 sample_group = c(rep("H1", 3), rep("H2",3), rep("H3",3), rep("H4",3), rep("H5",3), rep("H6",3)),
                 stringsAsFactors = FALSE)


Hseries <- readMsExperiment(spectraFiles = xcmsXML, sampleData = pd)

table(fromFile(Hseries))

basepeaks <- chromatogram(Hseries, aggregationFun = "max")

#Spectra Visuals#####

group_colours <- paste0(brewer.pal(6, "Set1") [1:6], "60") 
names(group_colours) <- c("H1", "H2", "H3", "H4", "H5", "H6")

plot(basepeaks, col = group_colours[sampleData(Hseries)$sample_group])

basepeaks_1 <- basepeaks[1, 1]
rtime(basepeaks_1) |> head()

intensity(basepeaks_1) |> head()

ioncount <- spectra(Hseries) |> 
  tic () |>
  split (f = fromFile (Hseries))
boxplot (ioncount, col = group_colours[sampleData(Hseries)$sample_group], ylab = "intensity", main = "Total Ion Count")

#Filter Rt from 0.5 s to 10 minutes (600s)
Hseries <- filterRt(Hseries, rt = c(0.5, 600))

#Make a subset to check the peaks before running centwave on entire set ####
#xchr <- findChromPeaks(Hseries, param = CentWaveParam(snthresh = 2))

#chromPeaks(xchr)
#sample_colors <- group_colours[xchr$sampleData]
#bg <- sample_colors[chromPeaks(xchr)[, "column"]]

#plot(xchr, col = group_colours, peakBg = bg)

#Run the centwave program ####

#If you know your exact peaks are clean, integrate = 2 , otherwise use integrate = 1
#outliers will create weird peak shapes for integrate = 2
cwp <- CentWaveParam(peakwidth = c(1, 100), noise = 1000,
                     prefilter = c(1, 2000), integrate = 1,
                     ppm = 10, 
                     mzdiff = 0.005,
                     snthresh = 3)

#Merge nearby peaks mainly if the peak is a long one and it gets cut off with the first filter
#minProp is default 0.75 and means that it will merge signals only if the large signal is higher than 75% of the small signal 
#mpp <- MergeNeighboringPeaksParam(expandRt = 3, minProp = 0.75)
#Hseries_refined <- refineChromPeaks(Hseries, mpp)

Hseries <- findChromPeaks(Hseries, param = cwp)

chromPeaks(Hseries) |>
  head()


#chromPeaks(xchr)

#neighbours <- MergeNeighboringPeaksParam (expandRt = 0.02)
#Hseries <- refineChromPeaks(Hseries, neighbours)


Hseries <- adjustRtime(Hseries, param = ObiwarpParam(binSize = 0.5))

pdp <- PeakDensityParam(sampleGroups = sampleData(Hseries)$sample_group,
                        minFraction = 0.5, minSamples = 2,
                        bw = 5, 
                        binSize = 0.015
)
Hseries <- groupChromPeaks(Hseries, param = pdp)

#pgp <- PeakGroupsParam(minFraction = 0.5, smooth = "loess", span = 0.4)

#Hseries <- adjustRtime(Hseries, param = pgp)


#Hseries <- groupChromPeaks(
#  dropFeatureDefinitions(Hseries),
#  PeakDensityParam(sampleGroups = sampleData(Hseries),
#                 minFraction = 0.5, bw = 5))

Hseries <- fillChromPeaks(Hseries, param = ChromPeakAreaParam())

featureValues(Hseries, value = "into") |> head()


grouped_peaks <- chromPeaks(Hseries)

grouped_peaks_df <- as.data.frame(grouped_peaks)

#write.csv(grouped_peaks_df, file = "Jourintegrate1_nostds_ppm2.csv", row.names = FALSE)


#library(SummarizedExperiment)

#Sum the integrated peaks -> can use maxint instead for the tallest integrated peak
res <- quantify(Hseries, value = "into", method = "sum")
rowData(res)
assay(res) |> head()

groupedassay <- as.data.frame(rowData(res))
groupedassay2 <- as.data.frame(assay(res))

write.csv(groupedassay, file = "groupedassayHseries1a.csv", row.names = TRUE)
write.csv(groupedassay2, file = "groupedassayHseries1b.csv", row.names = TRUE)


library(mzR)

writeMSData(
  Hseries,
  file,
  outformat = c("mzxml"),
  merge = FALSE,
  verbose = isMSnbaseVerbose(),
  copy = FALSE,
  software_processing = NULL
)

assays(res)$raw_nofill <- featureValues(Hseries, filled = FALSE, method = "sum")
assay(res, "raw_nofill") |> head()

write.csv(peaklist, "cactusPeaks2.csv")

#Visual exploration with plots and tables ####
library(tidyr)
library(dplyr)
library(readxl)
library(ggfortify)

cactusorg <- read_excel("xcmsHseries.xlsx")
cactusorg2 <- cactusorg[,c(1,7,10:27)]

peaklistlong <- pivot_longer(cactusorg2, cols=c("H1_1.mzXML","H1_2.mzXML", "H1_3.mzXML", "H2_1.mzXML","H2_2.mzXML", "H2_3.mzXML", "H3_1.mzXML","H3_2.mzXML", "H3_3.mzXML","H4_1.mzXML","H4_2.mzXML", "H4_3.mzXML","H5_1.mzXML","H5_2.mzXML", "H5_3.mzXML", "H6_1.mzXML","H6_2.mzXML", "H6_3.mzXML"), names_to = "sample", values_to = "peak_area")
str(peaklistlong)

sep1 <- strsplit ((as.character(peaklistlong$sample)),split = "_")
peaklistlong$sample <- sapply (sep1,"[[",1)

peaklistlong$sample <- as.factor(peaklistlong$sample)

#lowerRt<- peaklistlong[peaklistlong$mz > 212.128 & peaklistlong$mz < 212.1289,]
#lowerRt<- peaklistlong[peaklistlong$mz == 212.1280,]

pca_res <- prcomp(peaklistlong[c(1,2,4)], scale. = TRUE)

autoplot(pca_res)

autoplot(pca_res, data = peaklistlong, colour = "sample")

summary(pca_res)

var_pca <- pca_res$sdev^2
per_var_pca <- round((var_pca / sum(var_pca)) * 100, 2)

per_var_pca <- tibble(
  PC = 1:length(per_var_pca),
  PER_VAR = per_var_pca
)

bar_pca <- per_var_pca[1:3,] %>% 
  ggplot(aes(x = as.factor(PC), y = PER_VAR)) +
  geom_col() + 
  xlab("PC") +
  ylab("% of total variance")
bar_pca

par()              # view current settings
par(mfrow = c(1, 1), mar = c(5,5,3,3)) 


boxplot((log(peaklistlong$peak_area)) ~ peaklistlong$sample)

library(cluster)
autoplot(clara(peaklistlong[c(1,2,4)], 6))
autoplot(pam(peaklistlong[c(1,2,4)],6), data = peaklistlong, frame = TRUE, frame.type = 'norm', frame.colour = 'sample')

#https://stackoverflow.com/questions/75287022/change-legend-and-shape-in-ggbiplot-pca


ggplot(peaklistlong, aes(x = mzmed, y = peak_area, color = sample)) +
  geom_point()

ggplot(peaklistlong, aes(x = sample, y = peak_area, color = mzmed)) +
  geom_point()+
  theme_bw()

ggplot(peaklistlong, aes(x = sample, y = log(peak_area), color = mzmed)) +
  geom_point() +
  theme_bw()

hist(peaklistlong$peak_area)


cactus.lm = lm(logjour ~ logwill) 
cactus.res = resid(cactus.lm)


#Cactus Scatterplot ##### 
cactusscatter <- read_excel("xsmsHseries_H1H6.xlsx")

logH1 <- log(cactusscatter$H1_average)

cactusscatter$H1_average[cactusscatter$H1_average == 0] <- NA
logH6 <- log(cactusscatter$H6_average)

cactus.lm = lm(cactusscatter$H1_average ~ cactusscatter$H6_average) 
cactus.res = resid(cactus.lm)

plot(cactusscatter$H6_average, cactus.res, 
     ylab="Residuals", xlab="WillAve", 
     main="Residual Plot for Metabolite") 



qqnorm(logH1)
qqline(logH1, col = "red")

qqnorm(logH6)
qqline(logH6, col = "red")




plot(logH1,logH6,xlim = c(5,20),ylim = c(5,20),xlab = "H1 Log Peak Area",ylab = "H6 Log Peak Area")
arrows(-1000,-1000,10000000,10000000,col = "red")
arrows(log(0.4),-100,log(0.8),100,col = "orange")
arrows(-10,-0.8,10,0.8,col = "orange")

plot(cactusscatter$H1_average,cactusscatter$H6_average,xlim = c(-1,10000000),ylim = c(-1,10000000),xlab = "H1 Peak Area",ylab = "H6 Peak Area")
arrows(-1000,-1000,10000000,10000000,col = "red")


### Scatter 2 ####
cactusscatter2 <- read.csv("Metabolites/groupedassay3_3.csv", header = TRUE)

#cactusscatter2 <- cactusscatter[,c(-3:-14)]

#cactusscatter2$JourAve[is.na(cactusscatter2$JourAve)] <- 0.00001
#cactusscatter2$WillAve[is.na(cactusscatter2$WillAve)] <- 0.00001

cactusscatter2no_NA <- na.omit(cactusscatter2)

cactusscatter2no_NA$JourAve <- as.numeric(cactusscatter2no_NA$JourAve)
cactusscatter2no_NA$WillAve <- as.numeric(cactusscatter2no_NA$WillAve)

str(cactusscatter2no_NA)

annotationscatter <- cactusscatter2no_NA[cactusscatter2no_NA$ppm == 0.44080 | cactusscatter2no_NA$ppm == 0.91548 | cactusscatter2no_NA$ppm == 0.47296 | cactusscatter2no_NA$ppm == 0.67921 | cactusscatter2no_NA$ppm == 1.01962, ] 

str(cactusscatter2)
cactusscatter2no_NA[!is.numeric(cactusscatter2no_NA),]

logjour <- log(cactusscatter2no_NA$JourAve)
logwill <- log(cactusscatter2no_NA$WillAve)

cactus.lm = lm(logjour ~ logwill) 
cactus.res = resid(cactus.lm)

plot(logwill, cactus.res, 
     ylab="Residuals", xlab="WillAve", 
     main="Residual Plot for Metabolite") 



qqnorm(logjour)
qqline(logjour, col = "red")

qqnorm(logwill)
qqline(logwill, col = "red")

plot(cactusscatter2$Jour_Ave,cactusscatter2$Will_Ave,xlim = c(-1,25000000),ylim = c(-1,25000000),xlab = "Average Jourdaniana Peak Height",ylab = "Average Williamsii Peak Height",
     family="sans", 
     cex.lab = 1.2,
     cex.axis = 1.3,
     bty='n')
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)
arrows(-1000,-1000,100000000,100000000,col = "black")
arrows(10000,0,10000000,130000000,col = "black")
arrows(0,10000,130000000,15000000,col = "black")
points((annotationscatter$Jour_Ave), (annotationscatter$Will_Ave), pch = 19, col=col2)

legend(x = "topright",
       legend = matching2$annotation,
       col = col2,
       pch = 19,
       inset = -0.15)

library(dplyr)

matching <- read.table("~/Cactus/Metabolites/alkaloidsfromgroupedassay3.txt", header = T)

#Can't retain enough decimals for dplyr to recognize it
#alkaloidmatch <- cactusscatter2[cactusscatter2$JourAve %in% matching$jourave,]
#alkaloidmatch$enzyme <- matching$V2[match(bigmatch$V4, matching$V3)]
matching2 <- matching[c(2,6,10,15,16),]

matching2$logjour <- log(matching2$jourave)
matching2$logwill <- log(matching2$willave)

col2 <- c("red","yellow","green","purple", "blue")

plot(logjour,logwill,xlim = c(-1,20),ylim = c(-1,20),xlab = "Average Jourdaniana Log Peak Height",ylab = "Average Williamsii Log Peak Height",
     family="sans", 
     cex.lab = 1.2,
     cex.axis = 1.3,
     bty = 'n')
axis(side = 1, lwd = 1.2)
axis(side = 2, lwd = 1.2)
arrows(-1000,-1000,10000000,10000000,col = "black")
points((matching2$logjour), (matching2$logwill), pch = 19, col=col2)

#parorig <- par()
#par(mar = c(5,5,2,2))
plot(cactusscatter$Jour_Ave,cactusscatter$Will_Ave,xlim = c(-1,25000000),ylim = c(-1,25000000),xlab = "Average Jourdaniana Peak Height",ylab = "Average Williamsii Peak Height",
     family="sans", 
     cex.lab = 1.2,
     cex.axis = 1.3,
     bty='n')
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)
arrows(-1000,-1000,100000000,100000000,col = "black")
arrows(10000,0,10000000,130000000,col = "black")
arrows(0,10000,130000000,15000000,col = "black")
points((matching2$jourave), (matching2$willave), pch = 19, col=col2)
legend(x = "topright",
       legend = matching2$annotation,
       col = col2,
       pch = 19,
       inset = -0.15)
par (xpd = F)

plot(cactusscatter2no_NA$JourAve,cactusscatter2no_NA$WillAve,xlim = c(-1,25000000),ylim = c(-1,25000000),xlab = "Average Jourdaniana Peak Height",ylab = "Average Williamsii Peak Height",
     family="sans", 
     cex.lab = 1.2,
     cex.axis = 1.3,
     bty='n')
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)
arrows(-1000,-1000,100000000,100000000,col = "black")
arrows(10000,0,10000000,130000000,col = "black")
arrows(0,10000,130000000,15000000,col = "black")
points((matching2$jourave), (matching2$willave), pch = 19, col=col2)
legend(x = "topright",
       legend = matching2$annotation,
       col = col2,
       pch = 19,
       inset = -0.15)

plot(cactusscatter2no_NA$JourAve,cactusscatter2no_NA$WillAve,xlim = c(-1,100000),ylim = c(-1,100000),xlab = "Average Jourdaniana Peak Height",ylab = "Average Williamsii Peak Height")
arrows(-1000,-1000,10000000,10000000,col = "red")
arrows(0,0,100000,1300000,col = "orange")
arrows(0,0,1300000,150000,col = "orange")


matchingheat <- as.matrix(matching[, -1])
rownames(matchingheat) <- matchingheat[,1]


library ("gplots") 
library("RColorBrewer")
library("ggplot2")

heatmap.2 (matchingheat, 
           col = brewer.pal(n = 6, name = "Reds"), 
           scale = "column", 
           margins=c(8,8), 
           trace = "none", 
           cexCol = 1.1,
           cexRow = 1.1,
           Colv = NA, 
           Rowv = NA,
           key = T
           #sepwidth=c(0.05,0.1),
           #sepcolor = "white", #line colour it separates the columns with
           # rowsep=1:nrow(matchingheat),
           # colsep=1:ncol(matchingheat)
)

hist(cactusscatter2no_NA$JourAve, breaks = 5, main = "Jourdaniana Peak Height", xlab="Peak Heights", xlim=c(0,30000000))
hist(cactusscatter2no_NA$WillAve, breaks = 5, main = "Williamsii Peak Height", xlab="Peak Heights", xlim=c(0,30000000))
