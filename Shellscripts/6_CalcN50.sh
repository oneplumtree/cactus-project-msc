#!/bin/bash

#Trinity Salmon Abundance Estimation

# .sh file can be found in:
###################################################
#####Loop through directories to find fasta files
###################################################

for f1 in /data/poppy/domesticationpoppy/ncbihisat/zeno/*trinity*/;

do
        name=${f1%%-trinity.fasta} #only keep the front half of the name
        echo "Processing N50 stats for $name"
         /data/programs/trinityrnaseq-v2.14.0/util/TrinityStats.pl "$name"/Trinity-GG.fasta #> "$name"trinityN50stats.txt
done > zenotrinityN50Stats.txt #save them all in one file

