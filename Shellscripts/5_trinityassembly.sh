#!/bin/bash

##########################
####Trinity assembly
##########################

for f1 in /data/poppy/domesticationpoppy/trimgalore/zeno/*markdup.bam;

do
        name=${f1%%markdup.bam}
        echo "Processing $name file Trinity"
        /data/programs/trinityrnaseq-v2.14.0/Trinity --genome_guided_bam ${name}markdup.bam \
        --genome_guided_max_intron 10000 --max_memory 100G --CPU 10 \
        --output ${name}_trinity_out
done
