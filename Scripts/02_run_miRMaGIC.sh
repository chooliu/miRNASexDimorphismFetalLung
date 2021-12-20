# ==============================================================================
# 02_run_miRMaGIC.sh
# quantify .fastq --> miRNA counts; need miR-MaGiC, jre, snakemake
# ==============================================================================




# based on miRBase v22.1 reference files i'd previously generated
# see https://github.com/KechrisLab/miR-MaGiC

conda activate py3
cat > run_mirmagic.sh
#!/bin/bash
MIRNADIR=fetal_lung_miRNAseq
mkdir fetal_lung_miRNAseq/mirna_counts
for FOLDERNAME in `ls $MIRNADIR/trimmed_reads/`
do
  gunzip $MIRNADIR/trimmed_reads/$FOLDERNAME/*
  mkdir $MIRNADIR/mirna_counts/$FOLDERNAME
  for FILE in `ls $MIRNADIR/trimmed_reads/$FOLDERNAME/*.fastq`
    do 
    FILENAME=$(basename $MIRNADIR/trimmed_reads/$FOLDERNAME/$FILE)
    snakemake \
    --directory $HOME/miR-MaGiC/pipeline \
    --snakefile $HOME/miR-MaGiC/pipeline/Snakefile \
    --config \
    outdir=$MIRNADIR/mirna_counts/$FOLDERNAME/ \
    fastq=$MIRNADIR/trimmed_reads/$FOLDERNAME/$FILENAME \
    mirna=$HOME/miR-MaGiC/resources/reference_sequences/miRBase_v22.1_mature_sequences_Homo_sapiens.fasta \
    mirna_gp=$HOME/miR-MaGiC/resources/group_tables/miRBase_v22.1_group_by_core_ID_Homo_sapiens.txt \
    jar=$HOME/miR-MaGiC/pipeline \
    k=20 \
    plus_strand_only=True
    done
done

sh run_mirmagic.sh > run_mirmagic_out
