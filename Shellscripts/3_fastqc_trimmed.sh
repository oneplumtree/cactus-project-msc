#!/bin/bash

#########################
#####Fastqc Post-Trim
#########################

conda activate base

for f1 in /data/poppy/domesticationpoppy/trimgalore/otherpoppy/;
do
        dir=/data/poppy/domesticationpoppy/fastqc_trimmed/
        echo "Processing $f1 FastQC Post-Trim"
        /data/programs/FastQC/fastqc --format fastq -t 8 --outdir \
        "$dir" \
        "$f1"
        echo 
        echo 
done


#########################
#####MultiQC
#########################

conda activate multiqc

multiqc /data/poppy/domesticationpoppy/trimgalore/otherpoppy -o /data/poppy/domesticationpoppy/trimgalore/otherpoppy --filename poppydomestication_otherpoppytrimmed_fastqc.html --title "Other Poppy Trimmed FastQC"

#######Used WinSCP to transfer HTML files to computer and looked on local browser