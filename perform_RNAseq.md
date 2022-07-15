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

To do this, a complete reference file with fasta and gtf information.

You can get them from NCBI website (https://www.ncbi.nlm.nih.gov/refseq/)

```
STAR --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir /dartfs-hpc/rc/home/d/f0033vd/Star_index_USA300_FPR3757 \
  --genomeFastaFiles /dartfs-hpc/rc/home/d/f0033vd/USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic.fna \
  --sjdbGTFfile /dartfs-hpc/rc/home/d/f0033vd/USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic_with_sRNA.gtf \
  --genomeSAindexNbases 8
```

Once you make a reference, it is time to align your reads.
```
STAR --genomeDir /dartfs-hpc/rc/home/d/f0033vd/Star_index_USA300_FPR3757 \
  --readFilesIn /dartfs-hpc/scratch/cheung/JC/results/cutadapt/S1_R1.trim.fastq.gz /dartfs-hpc/scratch/cheung/JC/results/cutadapt/S1_R2.trim.fastq.gz \
  --readFilesCommand zcat \
  --sjdbGTFfile /dartfs-hpc/rc/home/d/f0033vd/USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic_with_sRNA.gtf \
  --runThreadN 4 \
  --alignIntronMax 1 \
  --outBAMsortingBinsN 200 \
  --limitBAMsortRAM 2000000000 \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterType BySJout \
  --outFileNamePrefix S1.
```

Bash file
```
#!/bin/bash
#SBATCH --job-name=multicore_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1:50:00
#SBATCH --mail-type=BEGIN,END,FAIL

for i in {13..24}; do \

STAR --genomeDir /dartfs-hpc/rc/home/d/f0033vd/Star_index_USA300_FPR3757 \
  --readFilesIn /dartfs-hpc/scratch/cheung/JC/results/trim/S${i}_R1.trim.fastq.gz /dartfs-hpc/scratch/cheung/JC/results/trim/S${i}_R2.trim.fastq.gz \
  --readFilesCommand zcat \
  --sjdbGTFfile /dartfs-hpc/rc/home/d/f0033vd/USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic_with_sRNA.gtf \
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
If you have multiple IDs with identical name in the file, you will see a wired result with uncounted features.



```
#!/bin/bash
#SBATCH --job-name=multicore_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:50:00
#SBATCH --mail-type=BEGIN,END,FAIL

featureCounts -T 4 -p -t CDS -a /dartfs-hpc/rc/home/d/f0033vd/USA300_FPR3757/GCF_000013465.1_ASM1346v1_genomic_with_sRNA.gtf -o /dartfs-hpc/scratch/cheung/JC/results/counts_USA300/totalcounts_AM.txt /dartfs-hpc/scratch/cheung/JC/results/alignment_USA300_FPR3757_AM/*.Aligned.sortedByCoord.out.bam
```




