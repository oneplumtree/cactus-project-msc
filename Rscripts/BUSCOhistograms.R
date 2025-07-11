#Facchini Lab Ginny - Jan 3 2024

rm(list=ls (all=T))
library (dplyr)
library (ggplot2)

Smallchromo <- read.table ("take3smallquerylargeDB.txt", header = F)
Largechromo <- read.table ("take3LargequerysmallDB.txt", header = F)

Smallchromo %>%
  group_by(V1) %>%
  summarise(max = max(V3, na.rm=TRUE)) %>%
  {. ->> Smallchromomax }

ggplot(Smallchromomax, aes(x=max)) + geom_histogram(bins=30, col=I("#F8766D")) + labs(x ='Percent Identity', y='Count', title='Qiushi Small Scaffold Query Super Scaffold') + theme_classic() + theme(plot.title = element_text(hjust = 0.5))

Largechromo %>%
  group_by(V1) %>%
  summarise(max = max(V3, na.rm=TRUE)) %>%
  {. ->> Largechromomax }

ggplot(Largechromomax, aes(x=max)) + geom_histogram(bins=30, col=I("#619CFF")) + labs(x ='Percent Identity', y='Count', title='Qiushi Super Scaffold Query Small Scaffold') + theme_classic() + theme(plot.title = element_text(hjust = 0.5))

Smallchromomax$Legend <- 'small scaffolds'
Largechromomax$Legend <- 'super scaffolds'
combolargesmall <- rbind(Smallchromomax, Largechromomax)
ggplot(combolargesmall, aes(x=max, color=Legend)) + geom_histogram(bins=30) + labs(x ='Percent Identity', y='Count') + theme_classic()

### Box plots
boxplot(max~Legend,
        data=combolargesmall,
        main="Busco Scaffolds Blast - Qiushi's Data",
        xlab="Scaffold Type",
        ylab="Percent Identity",
        col="orange",
        border="brown"
)

summary(Smallchromomax)
summary(Largechromomax)
