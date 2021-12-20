# ==============================================================================
# 16_final_prefiltering_and_dataexport.R
# filter duplicated samples -- based on scripts 05 EDA (esp. batch effects)
# && export data to be made public (Github repo)
# ==============================================================================




# final metadata tibble --------------------------------------------------------
# order by total counts then select first run_sample value

filter_subjects_nodupes <-
  phenodata_naive %>%
  bind_cols(seq = 1:nrow(.)) %>%
  arrange(-totcounts) %>%
  filter(!grepl("Rerun", run)) %>% # alt: filter by total counts only
  filter(!duplicated(Recruitment_ID)) %>%
  .$seq %>%
  `%in%`(1:nrow(phenodata_naive), .)

# results in final n = 298 independent samples
phenodata_final <-
  phenodata_naive %>%
  filter(filter_subjects_nodupes) # %>%

# re-calculate age factors to make sure still orthogonal
phenodata_final$AgeLin <- poly(phenodata_final$Age, degree = 2)[ , 1]
phenodata_final$AgeSq <- poly(phenodata_final$Age, degree = 2)[ , 2]




# final miRNA matrix -----------------------------------------------------------
# 901 --> 899 recalculating 25% rule

filter_mirnas_in_enough_subjects_nodupes <-
  counts_naive %>%
  select(-run_sample, -Recruitment_ID) %>%
  t() %>%
  .[ , filter_subjects_nodupes] %>%
  apply(., 1, function(x) { sum(x != 0) }) %>%
  `>=`(0.25*sum(filter_subjects_nodupes)) %>%
  names(.)[.]

counts_final <-
  counts_naive %>%
  select(all_of(filter_mirnas_in_enough_subjects_nodupes)) %>%
  t() %>%
  .[ , filter_subjects_nodupes] %>%
  set_colnames(phenodata_final$Recruitment_ID)




# annotations ------------------------------------------------------------------
# (for final results export & check chromosome) 

annotation_miRNAs <-
  tibble(miRNA = rownames(counts_final)) %>%
  left_join(duplicated_miRNA_names %>%
              transmute(miRNA, Sequence = sequence_dna,
                        AltName = alternative_names),
            by = "miRNA") %>%
  left_join(., mirna_chrom_locations, by = "miRNA")




# public data export to Github -------------------------------------------------

# counts
counts_final %>%
  as_tibble(rownames = "miRNA") %>%
  set_colnames(., c("miRNA", paste0("Sample", 1:nrow(.))) ) %>%
  write_tsv("./PublicData/miRNA_counts.tsv")

# metadata
phenodata_final %>%
  transmute(Sample = paste0("Sample", 1:nrow(.)),
            Age, AgeLin, AgeSq, SexMale = as.numeric(Sex),
            SmokeExposure = as.numeric(SmokeBinary),
            Run = as.factor(run) %>% 
              ( function(x) {
                levels(x) <- paste0("Run", 1:4)
                return(x)
              } )) %>%
  write_tsv("./PublicData/sample_metadata.tsv")



