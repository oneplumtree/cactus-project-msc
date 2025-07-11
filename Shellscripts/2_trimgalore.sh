!/bin/bash

#########################
#####TrimGalore
#########################

###if necessary make a new conda environment and install trim galore
###conda create -n cutadapt
###conda activate cutadapt
###conda install bioconda::trim-galore  ###program installation
###conda install conda-forge::pigz	###paralellization of runs

conda activate cutadapt 

for f1 in /data/poppy/domesticationpoppy/rawsequences/*_R1.fastq.gz;
do
        name=${f1%%_R1.fastq.gz} ###Keep the prefix
        echo "Processing $name file"
        trim_galore --phred33 --paired -j 6 --fastqc ${name}_R1.fastq.gz ${name}_R2.fastq.gz \ 
	--output_dir /data/poppy/domesticationpoppy/trimgalore/otherpoppy
done
