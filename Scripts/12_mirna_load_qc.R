# ==============================================================================
# 12_mirna_qc.R
# EDA on miRNAs from script 00f, combine samplelane values --> "alias IDs"
# ==============================================================================




# attach .fastq -> sample ID mapping -------------------------------------------

# mirmagic_count_table_mirnabysample = 2627 mirnas x 420 "run_sample" combos
# fastq_labels = 1482 .fastq files x 7 identifiers (names, run, sample ids)
# "run_sample" combo is distinct from biological "sample" (alias/recruitment ID)
# due to technical replicates of samples across runs

# set name of HBRNA technical controls (included run-to-run) all as "HBRNA"
fastq_labels <-
  fastq_labels %>%
  mutate(Recruitment_ID =
           if_else(grepl("HBR", Recruitment_ID), "HBRNA", Recruitment_ID)) 

# join counts with .fastq mapping data frame
mirmagic_counts_plus_ID <-
  mirmagic_count_table_samplebymirna %>%
  left_join(fastq_labels %>%
              select(run_sample, Recruitment_ID) %>%
              mutate(Recruitment_ID = if_else(grepl("HBR", Recruitment_ID),
                                              "HBRNA", Recruitment_ID)) %>%
              filter(!duplicated(run_sample)), by = "run_sample")

# check which samples associated with multiple run-by-sample .fastq files
samples_has_tech_reps <-
  mirmagic_counts_plus_ID$Recruitment_ID %>% table %>% sort %>% .[`!=`(., 1)]




# annotation of miRNAs ---------------------------------------------------------
# genomic position & sequence checking / filtering

# sequences of mature miRNAs from mirbase v22.1
mirna_sequences_hsa <-
  readLines("Data/mature.fa") %>%
  data.frame(
    orig_header = .[c(T, F)],
    sequence = .[c(F, T)]) %>%
  transmute(id = str_split(orig_header, " ") %>% map_chr(., ~ .[[1]]),
            accession = str_split(orig_header, " ") %>% map_chr(., ~ .[[2]]),
            sequence_dna = gsub("U", "T", sequence)) %>%
  filter(grepl("hsa", id)) %>%
  mutate(id = gsub(">", "", id))


# annotate miRNAs with chr location (& important for autosomal-only analyses)
mirna_chrom_locations <-
  fread("Data/hsa.gff3") %>%
  tidyr::separate(col = V9, sep = "Name\\=", into = c("discard", "tmp")) %>%
  separate(col = tmp, sep = ";", into = c("Name", "disc")) %>%
  transmute(miRNA = Name, Chr = `#`, Left = V4, Right = V5)

# 16 mature sequences shared by miRNAs
duplicated_miRNA_names <-
  mirna_sequences_hsa %>%
  arrange(sequence_dna, id) %>%
  group_by(sequence_dna) %>%
  summarize(miRNA = id %>% first,
            N_mirna = id %>% unique %>% length,
            alternative_names = unique(id) %>% .[-1] %>%
              paste0(., collapse = " / "))

# miR-maGiC will generate multiple entries for same sequence
# eg., hsa-miR-516b-3p and hsa-miR-516a-3p same mature sequence --
# same count vector appears 2x under 516a and 516b; take first miRNA ID if dupe
filter_mirnas_matureseqsame <-
  mirna_sequences_hsa %>%
  arrange(id) %>%
  filter(id %in% names(mirmagic_count_table_samplebymirna)) %>%
  filter(!duplicated(sequence_dna)) %>%
  .$id 

# omit duplicated miRNAs
mirmagic_count_table_samplebymirna <-
  mirmagic_count_table_samplebymirna %>%
  select(c("run_sample", all_of(filter_mirnas_matureseqsame)))




# misc miRNA QC ----------------------------------------------------------------
# additional filtering / descriptive statistics on miRNA prevalence

# for each miRNA, number of run_samples the miRNA is observed within
# parabola-shaped distribution, with most density at 0 to 10% of samples
# 2,078 / 2,585 mature miRNA sequences in db observed (i.e., non-zero) in data
# later filtering step: omit miRNAs observed in <20% of dataset
percent_subjects_obs <-
  mirmagic_count_table_samplebymirna %>%
  select(-1) %>%
  summarize_all(function(x) { sum(x != 0) / 419 }) %>%
  pivot_longer(data = ., cols = names(.)) %>%
  .$value * 100

sum(percent_subjects_obs != 0) # 2078 / 2584
length(percent_subjects_obs)

ggplot(data = NULL, aes(x = percent_subjects_obs)) +
  geom_histogram(color = "white", fill = "darkgrey") +
  xlab("% of Subjects in which miRNA Detected") +
  geom_text(aes(x = 25, y = 125), label = "omit miRNAs obs\nin <25% samples",
            color = "red", hjust = 0, size = 3.5) +
  theme_classic() +
  geom_vline(xintercept = 25, color = "red", lty = 2) +
  scale_y_continuous("Count of Mature miRNAs", expand = c(0, 0))


# number of miRNAs observed in each run_sample set
# median(num_mirnas_each_sample) = 785, left-tailed distribution
num_mirnas_each_sample <-
  mirmagic_count_table_mirnabysample %>%
  select(-1) %>%
  summarize_all(function(x) { sum(x != 0) }) %>%
  pivot_longer(data = ., cols = names(.)) %>%
  .$value

ggplot(data = NULL, aes(x = num_mirnas_each_sample)) +
  geom_histogram(color = "white", fill = "darkgrey") +
  theme_classic() +
  xlab("# of Unique miRNAs Observed in Each Run-by-Sample") +
  scale_y_continuous("Count of Subjects", expand = c(0, 0))

# total library counts by miRNA number detected
# not a clear 1:1 monotonic association with sequencing depth & miRNA detection
ggplot(data = NULL,
       aes(x = mirmagic_count_table_mirnabysample %>%
             select(-name) %>% colSums(),
           y = mirmagic_count_table_mirnabysample %>%
             select(-name) %>% `!=`(0) %>% colSums())) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 100000, lty = 2, color = "red") +
  geom_text(aes(x = 5e5, y = 50),
            label = "omit samples with < 1x10^5\n'mapped' miRNA counts",
            color = "red", hjust = 0, size = 3.5) +
  theme_classic() +
  scale_y_continuous("# Unique miRNAs Detected") +
  scale_x_continuous("Total Counts")




# check consistency of technical replicates ------------------------------------
# (note: some potentially more informative looks at tech reps in script 05)

# which samples have duplicates
mirmagic_counts_plus_ID$Recruitment_ID %>% table %>% .[`!=`(., 1)]

# among top miRNAs, check correlation from replicate to replicate
# select some limited amount to avoid cor inflation due to sparsity
num_mirnas_include <- 100

filter_most_common_miRNAs <-
  mirmagic_counts_plus_ID %>%
  select(-run_sample, -Recruitment_ID) %>%
  summarize_all(sum) %>% 
  pivot_longer(cols = 1:ncol(.)) %>%
  arrange(-value) %>%
  .$name %>%
  .[1:num_mirnas_include]

min_corrs_in_dupes <-
  sapply(names(samples_has_tech_reps),
       function(aliasid, mincountfilter = 100000) {

         Nreps_effective <-
           mirmagic_counts_plus_ID %>%
           filter(Recruitment_ID %in% aliasid) %>%
           filter(., rowSums(select(., -Recruitment_ID, -run_sample)) > mincountfilter) %>%
           nrow %>% print

         if (Nreps_effective > 1) {
         mirmagic_counts_plus_ID %>%
           filter(Recruitment_ID %in% aliasid) %>%
           filter(., rowSums(select(., -Recruitment_ID, -run_sample)) > mincountfilter) %>%
           select(filter_most_common_miRNAs) %>%
           t %>%
           cor(method = "spearman") %>%
           min(na.rm = T)
         } else {
           NA
         }
       }) %>%
  sort

# correlations pretty high between technical replicates,
# lowest correlation is Spearman r = 0.900 using top 100 miRNAs,
# HBRNA correlation 0.928
min_corrs_in_dupes %>% .[`<`(., 0.97)]
min_corrs_in_dupes %>% .[`<`(., 0.9)]

# for comparison, compare correlation between two random samples
# median corr ~0.86 with long left tail; all technical samples have higher cor
null_distr_of_corrs <-
  1:10000 %>%
  sapply(.,
         function(seed, mincountfilter = 100000) {
           set.seed(seed)
           mirmagic_counts_plus_ID %>%
             filter(Recruitment_ID %in% (
               sample(unique(mirmagic_counts_plus_ID$Recruitment_ID), 2))) %>%
             slice(., sample(1:nrow(.), size = nrow(.))) %>%
             filter(!duplicated(Recruitment_ID)) %>%
             filter(., rowSums(select(., -Recruitment_ID, -run_sample)) >
                      mincountfilter) %>%
             select(filter_most_common_miRNAs) %>%
             t %>%
             cor(method = "spearman") %>%
             min(na.rm = T)
         })

# show null versus observed correlations
ggplot(data = NULL,
       aes(x = min_corrs_in_dupes %>% .[setdiff(names(.), "HBRNA")])) +
  geom_histogram(color = "white", fill = "darkgrey") +
  geom_line(aes(x = null_distr_of_corrs), stat = "density", color = "red", lty = 2) +
  geom_text(aes(x = 0.7, y = 3), label = 'correlation between\nrandom samples',
            color = "red", size = 3.5, hjust = 1) +
  geom_point(aes(x = 0.94, y = 0.5), shape = 8, color = "blue") +
  geom_text(aes(x = 0.94, y = 2), label = "HBRNA", color = "blue") +
  theme_classic() +
  xlab("Across-Replicate Correlation, by Recruitment ID") +
  scale_y_continuous("Count of Alias IDs", expand = c(0, 0))
