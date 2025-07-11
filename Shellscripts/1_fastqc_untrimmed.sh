!/bin/bash

#########################
#####FastQC Pre-Trim
#########################

for f1 in /data/poppy/domesticationpoppy/*;
do
        dir=/data/poppy/domesticationpoppy/qc_untrimmed/
        echo "Processing $f1 file"
        /data/programs/FastQC/fastqc --format fastq --outdir \
        "$dir" \
        "$f1"
        echo 
        echo 
done




