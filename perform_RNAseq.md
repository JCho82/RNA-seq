# Perform RNA-seq with raw data

### Cut out adapters

'''cutadapt \
-o /dartfs-hpc/scratch/cheung/JC/results/trim/S1_R1.trim.fastq.gz \
-p /dartfs-hpc/scratch/cheung/JC/results/trim/S1_R2.trim.fastq.gz \
/dartfs-hpc/scratch/cheung/JC/raw/AM001_S1_R1_001.fastq.gz /dartfs-hpc/scratch/cheung/JC/raw/AM001_S1_R2_001.fastq.gz \
-m 1 --nextseq-trim=20 -j 4 > S1.cutadapt.report
'''
