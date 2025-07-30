#!/bin/bash

###################################################
#####Calculate the N50 of Assembled Transcriptomes
###################################################

for f1 in /data/poppy/domesticationpoppy/hisat2/otherpoppy/*/;

do
        name=${f1%%-trinity.fasta} #only keep the front half of the name
        echo "Processing N50 stats for $name"
         /data/programs/trinityrnaseq-v2.14.0/util/TrinityStats.pl "$name"/Trinity-GG.fasta > "$name"trinityN50stats.txt
done
