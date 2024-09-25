
rm(list=ls (all=T))
library (dplyr)
library (ggplot2)

##Williamsii
williamsiiSmall <- read.table ("~/R/Cactus/GlobalGenes/williamsiiSmallDB.txt", header = F)

williamsiiSmall %>%
  group_by(V1, V2, V4) %>%
  summarise(percent = max(V3, na.rm=TRUE)) %>%
  {. ->> williamsiiSmallQmax }

williamsiiSmallQmax %>%
  group_by(V1, V2, percent) %>%
  summarise(max = max(V4, na.rm=TRUE)) %>%
  {. ->> williamsiiSmallQmax2 }

williamsiiSmallunique <- williamsiiSmallQmax2[!duplicated(williamsiiSmallQmax2[,c('V1')]),]

summary(williamsiiSmallunique)

ggplot(williamsiiSmallunique, aes(x=percent)) + geom_histogram(bins=30, col=I("#F8766D")) + labs(x ='Percent Identity', y='Count', title='Williamsii Small Compartment') + theme_classic() + theme(plot.title = element_text(hjust = 0.5))


williamsiiLarge <- read.table ("williamsiiLargeDB", header = F)

williamsiiLarge %>%
  group_by(V1, V2, V4) %>%
  summarise(percent = max(V3, na.rm=TRUE)) %>%
  {. ->> williamsiiLargeQmax }

williamsiiLargeQmax %>%
  group_by(V1, V2, percent) %>%
  summarise(max = max(V4, na.rm=TRUE)) %>%
  {. ->> williamsiiLargeQmax2 }
williamsiiLargeunique <- williamsiiLargeQmax2[!duplicated(williamsiiLargeQmax2[,c('V1')]),]

#write.csv(williamsiiLargeunique,"~/R/Cactus/williamsiiLargeunique.csv")

#williamsiiLargeunique <- read.csv ("williamsiiLargeunique.csv", header = T)


ggplot(williamsiiLargeunique, aes(x=percent)) + geom_histogram(bins=30, col=I("#619CFF")) + labs(x ='Percent Identity', y='Count', title='Williamsii Large Compartment') + theme_classic() + theme(plot.title = element_text(hjust = 0.5))


williamsiiSmallunique$Legend <- 'small scaffolds'
williamsiiLargeunique$Legend <- 'super scaffolds'
williamsiicombo <- rbind(williamsiiSmallunique, williamsiiLargeunique)

summary(williamsiiSmallunique)
summary(williamsiiLargeunique)


#Jour
jourSmall <- read.table("Jour_smalldb.txt", header = F)

jourSmall %>%
  group_by(V1, V2, V4) %>%
  summarise(percent = max(V3, na.rm=TRUE)) %>%
  {. ->> jourSmallmax}

jourSmallmax %>%
  group_by(V1, V2, percent) %>%
  summarise(max = max(V4, na.rm=TRUE)) %>%
  {. ->> jourSmallmax2}


jourSmallunique <- jourSmallmax2[!duplicated(jourSmallmax2[,c('V1')]),]

ggplot(jourSmallunique, aes(x=percent)) + geom_histogram(bins=30, col=I("#F8766D")) + labs(x ='Percent Identity', y='Count', title='Jourdaniana Small Compartment') + theme_classic() + theme(plot.title = element_text(hjust = 0.5))

summary (jourSmallunique)


JourLarge <- read.table ("Jour_superdb.txt", header = F)

JourLarge %>%
  group_by(V1, V2, V4) %>%
  summarise(percent = max(V3, na.rm=TRUE)) %>%
  {. ->> JourLargeQmax }

JourLargeQmax %>%
  group_by(V1, V2, percent) %>%
  summarise(max = max(V4, na.rm=TRUE)) %>%
  {. ->> JourLargeQmax2 }

JourLargeunique <- JourLargeQmax2[!duplicated(JourLargeQmax2[,c('V1')]),]

summary(JourLargeunique)

ggplot(JourLargeunique, aes(x=percent)) + geom_histogram(bins=30, col=I("#619CFF")) + labs(x ='Percent Identity', y='Count', title='Jourdaniana Large Compartment') + theme_classic() + theme(plot.title = element_text(hjust = 0.5))


