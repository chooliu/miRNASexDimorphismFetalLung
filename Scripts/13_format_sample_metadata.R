# ==============================================================================
# 13_format_sample_metadata.R
# sample QC && use sample metadata in exprs set --> DESeq2 ready fmt
# ==============================================================================




# starting metadata tibble -----------------------------------------------------

# primary metadata file via microarray dataset
# (Bioconductor ExpressionSet format prepped by collaborators [A.T.K.?])
# 444 microarray samples x 50 variables, extract features also relevant to miRNA 

phenodata <-
  FetalLungeset.rma %>%
  phenoData() %>%
  .@data %>%
  as_tibble() %>%
  mutate(
    PostNatal = grepl("PN", `Pre.Age..Days.`),
    Age = `Pre.Age..Days.` %>% gsub("PN", "", .) %>%
      as.character() %>% as.numeric(),
    Age = if_else(PostNatal, Age + 9 * 30, Age),
    Pseudoglandular = Age < 17 * 7,
    Cotinine = `Cotinine..ng.g.placenta.`,
    SmokeBinary = `Cotinine..ng.g.placenta.` >= 7.5,
    Sex = Predict.Sex.from.PCA.of.Y.chrom.probes == "M",
    YearReceived = Date.Received %>% as.character() %>%
      str_split("/") %>% map_chr(~ .[3])
  )




# add extra information --------------------------------------------------------
# add genotype PCs via a forthcoming eQTL manuscript in preparation
# (paired genotype data, but not available on all samples)
# 444 x 56
phenodata %<>%
  left_join(
    x = .,
    y = read_tsv("Data/20191001_fetal_lung_genotype_PCs.tsv") %>%
      mutate(ID = as.character(ID)),
    by = c("ALIAS_ID" = "ID")
  )

# add RINs from miRNA preps - via Kansas City collaborators
# 444 x 57
phenodata %<>%
  left_join(
    x = .,
    y = read_csv("./Data/FetalLung_2013-05-30_478_sample info for microRNA_SELECTION_withextraction1.csv") %>%
      select(`Clin Pharm ID`, RIN) %>% filter(!duplicated(`Clin Pharm ID`)),
    by = c("ALIAS_ID" = "Clin Pharm ID")
  )




# filter pseudoglandular samples -----------------------------------------------
# 403 x 57

phenodata_filtered <-
  phenodata %>%
  filter(!duplicated(ALIAS_ID)) %>%  # microarray tech dupes 
  filter(!PostNatal) %>% # remove three post-natal
  filter(Pseudoglandular == 1) %>% # remove 21 outside of 7-17 week inclusive
  filter(!is.na(SmokeBinary)) %>% # missing placental cotinine
  filter(!ALIAS_ID %in%
           c("20042", "20631", "20225", "21985", "20735", "20774")) # [*]

# [*] paired 850K EPIC methylation dataset
# where these samples appeared to have different sex predictions,
# exclude just in case to err on mitigating sample mix-ups

# total of 22 samples with miRNA profiles excluded
setdiff(fastq_labels$Recruitment_ID,
        phenodata_filtered$ALIAS_ID) %>%
  sort %>%
  length

