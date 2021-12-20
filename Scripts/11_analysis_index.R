# ==============================================================================
# 11_analysis_index.R
# preliminaries, incl. importing libraries, external data sources
# also shows script workflow below
# ==============================================================================




# key libraries ----------------------------------------------------------------

library(tidyverse)
library(data.table)
library(magrittr)
library(writexl)
library(readxl)

library(cowplot)
library(DESeq2)
library(RUVSeq)
library(qvalue)
library(vegan)

library(viridis)
library(cowplot)
library(ggConvexHull)
library(ggbeeswarm)
library(ggridges)
library(ggthemes)
library(eulerr)
library(scales)

# also, the following two libraries called later (fxn name / scoping issues)
# library(miRNAtap) & library(biomaRt)




palette_color_sex <-
  c(Male = "maroon", yes = "maroon", `TRUE` = "maroon",
    Female = "black", `FALSE` = "black")

palette_shape_sex <-
  c(Male = 19, Male = 19, `TRUE` = 19,
    Female = 1, `FALSE` = 1)




# folder structure -------------------------------------------------------------

c("Scripts", "Data", "FiguresTables",
  "PublicData", "Results") %>% sapply(dir.create)




# data loading -----------------------------------------------------------------

# output of script 00f:
# contains miRNA counts & mapping between .fastq names --> alias_id
load("./Data/miRNA_counts_by_runsample.Rdata")

# from existing work - microarray data / sample metadata
load("./Data/FetalLung_5.30.2013_444samples_RMA_Pheno.rda")

# output date to prevent overwriting
currentdateiteration <- "20211109"





# analysis scripts in order ----------------------------------------------------

# Step 0) run scripts 01 to 03 first (bioinformatics: .fastq --> miR counts)

# Step 1) general dataset QC, cleaning, and modeling
# (applies to our work more broadly on this prenatal lung miRNA dataset)

source("./Scripts/12_mirna_load_qc.R") # miRNA quants --> sample IDs && QC/EDA
source("./Scripts/13_format_sample_metadata.R") # clean metadata for modeling
source("./Scripts/14_prelim_sample_QC_naive.R") # prep samples for prelim EDA
source("./Scripts/15_deseq_models_naive.R") # DESeq2 w/ dupes, no RUV
source("./Scripts/16_final_prefiltering_and_dataexport.R") # final QC, PublicData export [***]
source("./Scripts/17_mod1.R") # RUVSeq & main effects model
source("./Scripts/18_mod2.R") # add age-by-sex interaction

# Step 2) analyses specific to sex analyses in this manuscript
source("./Scripts/21_sexeffect_primary_model.R") # diff exp miRNA by sex
source("./Scripts/22_agesex_interaction_effect.R") # sex-specific trajectories
source("./Scripts/23_interact_model_sensanalysis.R") # sens analysis for 11.R
source("./Scripts/24_make_figures.R") # figures in main text
source("./Scripts/25_run_and_format_miEAA.R") # export for GSEA 
source("./Scripts/26_target_annotation_and_export.R") # predicted targets, export final results tables


# [***] note on /PublicData/ folder:
# unfortunately, can't make all files needed to fully run scripts
# 12 through 16 public at this time -- however, the metadata data frame &
# miRNA counts matrix from this step can be used to fully reproduce 17 to 26
