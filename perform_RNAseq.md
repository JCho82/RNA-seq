# Perform RNA-seq with raw data

### Cut out adapters

```cutadapt \
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



