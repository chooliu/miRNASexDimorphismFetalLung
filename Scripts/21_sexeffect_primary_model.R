# ==============================================================================
# 21_sexeffect_primary_model.R 
# differences in mean miRNA levels between male and female samples
# format/process results from from primary model (mod1: script 17)
# ==============================================================================




# re-calculate DE by sex, after removing sex chromosomes -----------------------
# focus on miRNAs mapping to autosomal chromosomes -- not removed further
#  upstream b/c (i) presumed may be useful for DESeq2 (e.g., lib size; mean-var)
# & (ii) want export as part of core prenatal lung miRNA dataset

DEresults_sex <-
  resultsmod1_sex %>%
  results(tidy = T) %>%
  arrange(pvalue) %>%
  filter(!is.na(padj)) %>%
  left_join(annotation_miRNAs %>% filter(!duplicated(miRNA)),
            by = c("row" = "miRNA"))

DEresults_sex <-
  DEresults_sex %>%
  filter(!Chr %in% c("chrX", "chrY"))

DEresults_sex$qval <-
  qvalue(DEresults_sex$pvalue)$qvalue




# format for publication -------------------------------------------------------

FinalResults_Sex <-
  DEresults_sex %>%
  arrange(pvalue) %>%
  rowwise() %>%
  transmute(FC = 2^log2FoldChange,
            `Higher Levels In` =
              if_else(log2FoldChange > 0, "Males", "Females"),
            miRNA = row,
            `Effect (log2FC)` = log2FoldChange,
            `SE (log2FC)` = lfcSE,
            `LRT Statistic` = stat,
            `p-value` = pvalue,
            `q-value` = qval,
            Signif = if_else(`q-value` < 0.05, "*", "")) %>%
  mutate_if(is.numeric, format, digits = 2) 

Table_S1 <-
  FinalResults_Sex %>%
  filter(Signif == "*")
