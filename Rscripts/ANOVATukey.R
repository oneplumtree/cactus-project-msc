rm(list=ls(all=TRUE))
#install.packages("multcomp")
library(multcomp)

orbitrap <-read.table("H1_C1_orbitrapmasses.txt", header = T)
orbitrap_long <- gather(orbitrap, Standard, Micrograms, ThreeFourDM5HPEA:Anhalinine, factor_key=TRUE)


for(i in unique(orbitrap_long$Standard)){
  print(i)
  subsetsample <- orbitrap_long[orbitrap_long$Standard == i,]

#Check if the data is normally distributed can log transform them if not
  plot(qqnorm(subsetsample$Micrograms, pch = 1, frame = FALSE))
  plot(qqline(subsetsample$Micrograms, col = "steelblue", lwd = 2))
  print(shapiro.test(subsetsample$Micrograms))

 #Visual Analysis of the quartiles
  boxplot(Micrograms ~ Sample, data = subsetsample)
  str(subsetsample$Sample)
  subsetsample$Sample <- as.factor(subsetsample$Sample)
#ANOVA linear model from baseR 
  aovmodel <- aov(Micrograms ~ Sample, subsetsample)
  print(summary(aovmodel))
 #Tukey analysis based on the previously generated ANOVA
  turkey<- glht(aovmodel, linfct = mcp (Sample = "Tukey"))
  cld(turkey)
  print(cld(turkey))
  plot(TukeyHSD(aovmodel, conf.level = 0.95))
}
