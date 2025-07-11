#!/bin/bash

#########################
#####hisat2
#########################

##### Note the reference was already built for a different project, make sure to check if there is reference when starting

##### Hisat2 genome mapping
for f1 in /data/poppy/domesticationpoppy/trimgalore/otherpoppy/*R1_val_1.fq.gz;
do
        name=${f1%%R1_val_1.fq.gz}
        short=${name#*_i5.} #only keeps the name after "_i5."
        short="${short%_}"  #removes the hanging underscore at the end
        echo $short
        echo "Processing $name file with Hisat2"
        hisat2 -x /data/poppy/ginny/hisat2/poppygenomeHisat2 -p 8 \
	--summary-file -1 "$name"R1_val_1.fq.gz -2 "$name"R2_val_2.fq.gz
done

##########################
########sort the hisat files with samtools
##########################

for f1 in /data/poppy/domesticationpoppy/trimgalore/otherpoppy/*.hisat;
do
        name=${f1%%.hisat}
        samtools view -bS "$name".hisat > "$name".bam                                                                           	samtools sort -@ 8 "$name".bam -o "$name".sorted.bam
        samtools rmdup -S "$name".sorted.bam "$name".sorted.markdup.bam
        stringtie -o "$name".gtf "$name".sorted.markdup.bam
done

##########################



