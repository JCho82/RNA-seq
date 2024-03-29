# Perform RNA-seq with raw data

### Cut out adapters

```
cutadapt \
-o /results/trim/S1_R1.trim.fastq.gz \
-p /results/trim/S1_R2.trim.fastq.gz \
/raw/AM001_S1_R1_001.fastq.gz /raw/AM001_S1_R2_001.fastq.gz \
-m 1 --nextseq-trim=20 -j 4 > S1.cutadapt.report
```
The option "-m 1" indicates that you want to set the minimum length of reads to 1, meaning zero-length sequences are disregarded.

"--nextseq-trim=20" is uniquely used for illumina nextseq or Novaseq. 

If you are using other instruments, you can set the quality cut with option -q (for instance, -q 20).

-o and -p indicate output files.

Following them, you should write your target raw files. In this examples, two files (R1 and R2) are provided (paired-end reads).

You can make a bash file for this.

```
#!/bin/bash
#SBATCH --job-name=multicore_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL

for i in {1..10}; do \

cutadapt \
-o /results/trim/S${i}_R1.trim.fastq.gz \
-p /results/trim/S${i}_R2.trim.fastq.gz \
/raw/JC0${i}_S${i}_R1.fastq.gz /raw/JC0${i}_S${i}_R2.fastq.gz \
-m 1 --nextseq-trim=20 -j 4 > /results/trim/S${i}.cutadapt.report \

echo $i is completed for trimming.

done
```
The example above is showing a run with 10 samples named JC01_S1_R1.fastq.gz and JC01_S1_R2.fastq.gz with increment numbers. 

The output files will be S1_R1.trim.fastq.gz, S1_R2.trim.fastq.gz, and etc.

### Align with STAR

First, you need to make a reference with STAR.

To do this, find a complete reference file with fasta and gtf information.

You can get them from NCBI website (https://www.ncbi.nlm.nih.gov/refseq/)

```
mkdir /reference/Star_index_USA300_FPR3757
```
And then you need to make a folder to save references generated by STAR.

```
STAR --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir /reference/Star_index_USA300_FPR3757 \
  --genomeFastaFiles /USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic.fna \
  --sjdbGTFfile /USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic_with_sRNA.gtf \
  --genomeSAindexNbases 8
```
Everything is straightforward except genomeSAindexNbases.

The number can be calculated by {log2(genome size)/2} - 1.

The default is set to 14 but it would be better to check how people are setting for this.

Staphy aureus has 2.8M genome size, generating 9 but people are using 8 with s. aureus (look at the link below).

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4452430

Once you make a reference, it is time to align your reads.
```
STAR --genomeDir /reference/Star_index_USA300_FPR3757 \
  --readFilesIn /results/trim/S1_R1.trim.fastq.gz /results/trim/S1_R2.trim.fastq.gz \
  --readFilesCommand zcat \
  --sjdbGTFfile /USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic_with_sRNA.gtf \
  --runThreadN 4 \
  --alignIntronMax 1 \
  --outBAMsortingBinsN 200 \
  --limitBAMsortRAM 2000000000 \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterType BySJout \
  --outFileNamePrefix S1.
```
Since I am working on bacteria, I set alignIntronMaz to 1 (bacteria don't have introns).

For limitBAMsortRAM, I set it to 2000000000 since this amount of memory was allowed me to use.

The option "outBAMsortingBinsN" indicates the number of genome bins fo coordinate-sorting.
To reduce the memory stress in the shared computer, I set it to 200 (the default is 50).


For Bash file
```
#!/bin/bash
#SBATCH --job-name=multicore_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1:50:00
#SBATCH --mail-type=BEGIN,END,FAIL

for i in {1..10}; do \

STAR --genomeDir /reference/Star_index_USA300_FPR3757 \
  --readFilesIn /results/trim/S${i}_R1.trim.fastq.gz /results/trim/S${i}_R2.trim.fastq.gz \
  --readFilesCommand zcat \
  --sjdbGTFfile /USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic_with_sRNA.gtf \
  --runThreadN 4 \
  --alignIntronMax 1 \
  --outBAMsortingBinsN 200 \
  --limitBAMsortRAM 2000000000 \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterType BySJout \
  --outFileNamePrefix S${i}. \

echo sample$i is completed for alignment

done

```

### Count reads in features with FeatureCounts

You should carefully check if duplicates or triplicates of ID in gtf file are present.
If you have multiple IDs with the identical name in the file, you will see a wired result with uncounted features.



```
#!/bin/bash
#SBATCH --job-name=multicore_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:50:00
#SBATCH --mail-type=BEGIN,END,FAIL

featureCounts -T 4 -p -t CDS -a /USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic_with_sRNA.gtf -o /results/totalcounts.txt /results/*.Aligned.sortedByCoord.out.bam
```
-T (capitalized) indicates the number of CPU threads.

-p is used for paired-end reads. 

-t is to specify the feature type (GTF.featureType)

-g is to specify the attribute type used to group features (eg. exons) into meta-features (eg. genes) when GTF annotation is provided. ‘gene id’ by default.

*Since -g is not correctly working, you need to use -t option if need to change the feature type.






