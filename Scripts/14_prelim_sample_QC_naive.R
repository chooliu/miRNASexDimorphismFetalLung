# ==============================================================================
# 14_naive_sample_QC_naive.R
# some naive preliminary filtering for a "naive" analysis for EDA
# (i.e., modeling with no RUVSeq components, does not omit technical replicates) 
# ==============================================================================




# sample + miRNA filtering -----------------------------------------------------

# mapped miRNA library size (total counts)
# 419 --> 402 compiled .fastq with >~ 5x10^5
filter_subjects_totcount <-
  mirmagic_counts_plus_ID %>%
  select(-run_sample, -Recruitment_ID) %>%
  rowSums() %>%
  `>=`(400000) # adds one sample with 499,412

# 419 --> 390 has sample metadata
filter_subjects_haspheno <-
  mirmagic_counts_plus_ID$Recruitment_ID %in%
  phenodata_filtered$ALIASID

# 2586 --> 930 miRNAs in 25% of participants
filter_mirnas_in_enough_subjects <-
  mirmagic_counts_plus_ID %>%
  select(-run_sample, -Recruitment_ID) %>%
  apply(., 2, function(x) { sum(x != 0) }) %>%
  `>=`(0.25*nrow(mirmagic_counts_plus_ID)) %>%
  names(.)[.]




# apply filtering --------------------------------------------------------------
# 375 x 901
counts_naive <-
  mirmagic_counts_plus_ID %>%
  filter(filter_subjects_totcount & filter_subjects_haspheno) %>%
  select(c("run_sample",
           all_of(intersect(filter_mirnas_in_enough_subjects,
                            filter_mirnas_matureseqsame)),
           "Recruitment_ID"))

# get phenodata, make sure in same order as count dat
phenodata_naive <-
  counts_naive %>%
  select(run_sample, Recruitment_ID) %>%
  left_join(., fastq_labels %>% select(run_sample, run) %>%
              filter(!duplicated(run_sample))) %>%
  left_join(., phenodata_filtered %>% filter(!duplicated(ALIAS_ID)),
            by = c("Recruitment_ID" = "ALIAS_ID"))

phenodata_naive$AgeLin <- poly(phenodata_naive$Age, degree = 2)[ , 1]
phenodata_naive$AgeSq <- poly(phenodata_naive$Age, degree = 2)[ , 2]




# check filtering worked, order samples ----------------------------------------

dim(phenodata_naive)
identical(counts_naive$Recruitment_ID, phenodata_naive$Recruitment_ID)
order_by_recruitment <-
  order(counts_naive$Recruitment_ID)

counts_naive <- counts_naive[order_by_recruitment, ]
phenodata_naive <- phenodata_naive[order_by_recruitment, ]

phenodata_naive$totcounts <-
  counts_naive %>% select(-run_sample, -Recruitment_ID) %>% apply(., 1, sum)

