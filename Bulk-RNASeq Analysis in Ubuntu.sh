#!/bin/bash
~/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-3 SRRXXXXXXXX.sra -O.
fastqc -f fastq SRRXXXXXXXX_1.fastq SRRXXXXXXXX_2.fastq
mv SRRXXXXXXXX_1.fastq ~/ngs/samples
mv SRRXXXXXXXX_2.fastq ~/ngs/samples
cd ~/ngs/samples/
java -jar /usr/share/java/trimmomatic-0.39.jar PE SRRXXXXXXXX_1.fastq SRRXXXXXXXX_2.fastq ../trim/SRRXXXXXXXX_trimmed_1.fastq ../trim/SRRXXXXXXXX_unpaired_1.fastq ../trim/SRRXXXXXXXX_trimmed_2.fastq ../trim/SRRXXXXXXXX_unpaired_2.fastq LEADING:15 TRAILING:15 SLIDINGWINDOW:4:25 MINLEN:50
cd ~/ngs/hisat/
~/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1/hisat2 -q -x ../indexes/grch38_genome/grch38/genome -1 ../trim/SRRXXXXXXXX_trimmed_1.fastq -2 ../trim/SRRXXXXXXXX_trimmed_2.fastq --add-chrname -S SRRXXXXXXXX.sam
htseq-count SRRXXXXXXXX.sam hg38.ncbiRefSeq.gtf > SRRXXXXXXXX.count
cut -f 2 SRRXXXXXXXX.count | grep -Evc '^0$'
cd


