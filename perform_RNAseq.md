# Perform RNA-seq with raw data

### Cut out adapters

```
cutadapt \
-o /dartfs-hpc/scratch/cheung/JC/results/trim/S1_R1.trim.fastq.gz \
-p /dartfs-hpc/scratch/cheung/JC/results/trim/S1_R2.trim.fastq.gz \
/dartfs-hpc/scratch/cheung/JC/raw/AM001_S1_R1_001.fastq.gz /dartfs-hpc/scratch/cheung/JC/raw/AM001_S1_R2_001.fastq.gz \
-m 1 --nextseq-trim=20 -j 4 > S1.cutadapt.report
```
The option "-m 1" is set to the minimum length of reads as 1
"--nextseq-trim=20" is uniquely used for illumina nextseq or Novaseq. 
If you are using other instruments, you can set the quality cut as "-q 20" 

-o and -p optinos for output files.

You can make bash file for this.

```
#!/bin/bash
#SBATCH --job-name=multicore_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL

for i in {13..24}; do \

cutadapt \
-o /dartfs-hpc/scratch/cheung/JC/results/trim/S${i}_R1.trim.fastq.gz \
-p /dartfs-hpc/scratch/cheung/JC/results/trim/S${i}_R2.trim.fastq.gz \
/dartfs-hpc/scratch/cheung/JC/raw/AM0${i}_S${i}_R1_001.fastq.gz /dartfs-hpc/scratch/cheung/JC/raw/AM0${i}_S${i}_R2_001.fastq.gz \
-m 1 --nextseq-trim=20 -j 4 > /dartfs-hpc/scratch/cheung/JC/results/trim/S${i}.cutadapt.report \

echo $i is completed for trimming.

done
```

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

Once you make the reference, it is time to align your reads.
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

