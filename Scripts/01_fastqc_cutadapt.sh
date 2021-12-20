# ==============================================================================
# 01_fastqc_cutadapt.sh
# check .fastq quality / presence of artifactual sequences --> trim + QC
# run on Centos6; prereqs include conda, py3, R, FastQC
# ==============================================================================




# run FastQC on raw .fastq files from sequencing core --------------------------
# .fastqs under ./fetal_lung_miRNAseq/raw_reads/[run name]

cat > runfastqc_mirna_raw.sh
#!/bin/bash
FASTQC=$HOME/FastQC/fastqc
TMPFOLDER=tmp
MIRNADIR=fetal_lung_miRNAseq
mkdir $MIRNADIR/fastqc_raw

for FOLDER in `ls $MIRNADIR/raw_reads/fastq`
do 
mkdir $MIRNADIR/fastqc_raw/$FOLDER/ 
  $FASTQC `ls $MIRNADIR/raw_reads/fastq/$FOLDER/*.fastq.gz` -o \
  $MIRNADIR/fastqc_raw/$FOLDER/ --noextract -d $TMPFOLDER -t 12
done




# generate cutadapt commands ---------------------------------------------------
# (used R to loop through all files due to some particulars of our server setup)

R
library(tidyverse)
setwd('fetal_lung_miRNAseq')
system('conda activate py3')
system('find raw_reads -name "*.fastq.gz" > raw_fastq_file_locations')

rawdat_filepath <- read_lines("raw_fastq_file_locations")
run_folder <-
  dirname(rawdat_filepath) %>%
  gsub("raw_reads\\/fastq\\/|raw_reads\\/|\\/Files", "", .)
run_folder %>% unique %>%
  paste0("fetal_lung_miRNAseq/trimmed_reads/", .) %>% sapply(dir.create)
outname <-
  paste0("fetal_lung_miRNAseq/trimmed_reads/",
         run_folder, "/", basename(rawdat_filepath))
outname %>% .[duplicated(.)] %>% any

cutadaptcmds <-
  tibble(rawdat_filepath,
         run_folder,
         outname) %>%
  mutate(cmd =
           paste("~/.local/bin/cutadapt -u 4 -a NNNNTGGAATTCTCGGGTGCCAAGG -j 48 -q 20 -m 10 --trim-n -o",
                 outname, rawdat_filepath, "; sleep 15;"))

write_lines(cutadaptcmds$cmd,
            path = "Scripts/run_cutadapt.sh")

system("sh Scripts/run_cutadapt.sh > Scripts/run_cutadapt.log")




# run FastQC on trimmed files --------------------------------------------------

cat > runfastqc_mirna_trimmed.sh
#!/bin/bash
FASTQC=$HOME/FastQC/fastqc
TMPFOLDER=tmp
MIRNADIR=fetal_lung_miRNAseq
mkdir $MIRNADIR/fastqc_trimmed

for FOLDER in `ls $MIRNADIR/trimmed_reads/`
do 
mkdir $MIRNADIR/fastqc_trimmed/$FOLDER/ 
  $FASTQC `ls $MIRNADIR/trimmed_reads/$FOLDER/*.fastq.gz` -o \
  $MIRNADIR/fastqc_trimmed/$FOLDER/ --noextract -d $TMPFOLDER -t 12
done
