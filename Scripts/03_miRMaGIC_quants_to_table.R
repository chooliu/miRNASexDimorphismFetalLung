# ==============================================================================
# 03_miRMaGIC_quants_to_table.R
# parse output of script 02 into count table
# ==============================================================================




# run in R  --------------------------------------------------------------------
setwd('fetal_lung_miRNAseq/mirna_counts')
library(tidyverse)
library(readxl)




# get filepaths & examine file format ------------------------------------------

mirmagic_output_paths <-
  list.files(recursive = T, full.names = T, pattern = ".txt")

# check format of files
# confirm that miR-MaGIC outputs the miRNA IDs all in same order
mirmagic_count_table <-
  map(mirmagic_output_paths, function(x) {
    read_tsv(x, col_names = F) %>% set_names(., c("mirna", x))
    })
identical(mirmagic_count_table[[1]]$mirna,
          mirmagic_count_table[[1000]]$mirna)


# because miRNAs in same order, load mirna ids first, then counts separately
mirna_ids <-
  read_tsv(mirmagic_output_paths[1], col_names = F, col_types = "c_") %>%
  set_names(., "mirna_name")

# should be 2627 x 1483
mirmagic_count_table <-
  map(mirmagic_output_paths, function(x) {
    read_tsv(x, col_names = F, col_types = "_i") %>%
      set_names(., gsub("\\.\\/|final_counts\\.|\\.txt", "", x))
  }) %>%
  bind_cols(mirna_ids, .)

# check format
mirmagic_count_table[1:4, 1:4]




# sample info  -----------------------------------------------------------------
fastq_labels <-
  enframe(mirmagic_count_table %>% names %>% .[-1],
          name = "index", value = "run_sample_lane") %>%
  mutate(run = run_sample_lane %>% str_split(., pattern = "/") %>% map_chr( ~ .[1]),
         sample_lane =  run_sample_lane %>% str_split(., pattern = "/") %>% map_chr( ~ .[2]),
         sample = sample_lane %>% str_split(., pattern = "_L") %>% map_chr( ~ .[1]),
         run_sample = paste0(run, "__", sample))

fastq_labels %>%
  summarize_all(function(x) { duplicated(x) %>% any })




# load sample manifest, maps .fastq --> ALIAS_ID -------------------------------
# thanks to ATK & RC (collaborators at Channing)

channing_manifest <-
  read_excel("../metadata/FetalLung_2019-09-19_QC.xls")

# for .fastq --> ALIAS_ID
channing_manifest <-
  channing_manifest %>%
  mutate(run_sample_lane =
           paste0(directory, "/",
                  gsub("\\.fastq\\.gz", "", filename_extract_from_fullpath))) 

# only samples not matched by "run_sample_ are HBRNA (technical reps)
# (in manifest, filepath is "[run]/Samples/HBRNA/" with run context;
# whereas mine are stored in "/HBRNA_fastq/ without run information")
channing_manifest$run_sample_lane %>% .[!`%in%`(., fastq_labels$run_sample_lane)]

channing_manifest <-
  channing_manifest %>%
  mutate(run_sample_lane_edit = 
           if_else(grepl("HBRNA", run_sample_lane),
                   paste0("HBRNA_fastq/",
                          gsub("\\.fastq\\.gz", "", filename_extract_from_fullpath)),
                   run_sample_lane))

# this edit to HBRNA file paths then make everything joinable on run_sample_lane
channing_manifest$run_sample_lane_edit  %>% .[!`%in%`(., fastq_labels$run_sample_lane)]

fastq_labels <-
  fastq_labels %>%
  left_join(channing_manifest %>% select(Recruitment_ID, run_sample_lane_edit),
            by = c("run_sample_lane" = "run_sample_lane_edit"))




# generate count table ---------------------------------------------------------

# 419 unique run x sample combinations
fastq_labels %>%
  group_by(run_sample) %>%
  group_keys %>%
  nrow

mirmagic_count_table_transpose <-
  pivot_longer(mirmagic_count_table,
               cols = names(mirmagic_count_table) %>% setdiff("mirna_name")) %>%
  pivot_wider(id_cols = "name", names_from = "mirna_name")

names(mirmagic_count_table_transpose)[1] <- "run_sample_lane"

# check some random entires to see if correct : should return 442
mirmagic_count_table %>%
  filter(mirna_name == "hsa-miR-99b-3p") %>%
  select(`10537388472-0FetalLungPT2/S-001614520_S2_L003_R1_001`)

mirmagic_count_table_transpose %>%
  filter(run_sample_lane == "10537388472-0FetalLungPT2/S-001614520_S2_L003_R1_001") %>%
  select(`hsa-miR-99b-3p`)


# finally, join the count matrix with fastq_labels & aggregate by run_sample (exclude lane)
# reasonable because no apparent lane-to-lane effect (previously examined; code not shown)
# 419 x 2,628
mirmagic_count_table_samplebymirna <-
  mirmagic_count_table_transpose %>%
  left_join(., fastq_labels %>%
              select(run_sample, run_sample_lane),
            by = "run_sample_lane") %>%
  group_by(run_sample) %>%
  summarize_at(-1, sum)  # ** can replace this with metrics to check lane-to-lane variance

# 2,628 x 419
mirmagic_count_table_mirnabysample <-
  pivot_longer(mirmagic_count_table_samplebymirna,
               cols = names(mirmagic_count_table_samplebymirna) %>% setdiff("run_sample")) %>%
  pivot_wider(id_cols = "name", names_from = "run_sample")




# final export is miRNA count table, aggregated by sample_run ------------------
# input to script 11
save(mirmagic_count_table_samplebymirna,
     mirmagic_count_table_mirnabysample,
     fastq_labels,
     file = "miRNA_counts_by_runsample.Rdata")


