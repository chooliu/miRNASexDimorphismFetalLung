# ==============================================================================
# 22_agesex_interaction_effect.R
# differences in miRNA age-trajectories between male and female samples
# format/process results from age-by-sex interaction model (mod2: script 18)
# ==============================================================================




# re-calculate DE after removing sex chromosomes -------------------------------
# interaction effect

DEresults_interaction <-
  resultsmod2_firstorderinteract %>%
  results(tidy = T) %>%
  arrange(pvalue) %>%
  filter(!is.na(padj)) %>%
  left_join(annotation_miRNAs %>% filter(!duplicated(miRNA)),
            by = c("row" = "miRNA"))

DEresults_interaction <-
  DEresults_interaction %>%
  filter(!Chr %in% c("chrX", "chrY"))

DEresults_interaction$qval <-
  qvalue(DEresults_interaction$pvalue)$qvalue




# format for publication -------------------------------------------------------

FinalResults_interaction <-
  DEresults_interaction %>%
  arrange(pvalue) %>%
  transmute(miRNA = row,
            `Slope Higher In` = if_else(log2FoldChange > 0, "Males", "Females"),
            `Effect (log2FC)` = log2FoldChange,
            `SE (log2FC)` = lfcSE,
            `LRT Statistic` = stat,
            `p-value` = pvalue,
            `q-value` = qval,
            Signif = if_else(`q-value` < 0.05, "*", "")) %>%
  mutate_if(is.numeric, format, digits = 2) 


Table_S3 <-
  FinalResults_interaction %>%
  filter(Signif == "*")




# check specifics of male-only and female-only age trajectories ----------------
# test if slope is signif. non-zero and slope direction (incr, decr)
# add the directionality of each to the final table

# female samples

DEresults_female <-
  resultsmod2_femalelinest %>%
  arrange(pvalue) %>%
  filter(row %in% DEresults_interaction$row) %>%
  left_join(annotation_miRNAs, by = c("row" = "miRNA"))

DEresults_female <-
  DEresults_female %>%
  filter(!Chr %in% c("chrX", "chrY"))

DEresults_female$qval <-
  qvalue(DEresults_female$pvalue)$qvalue

FinalResults_female <-
  DEresults_female %>%
  arrange(pvalue) %>%
  transmute(miRNA = row,
            `Trajectory` = if_else(log2FoldChange > 0, "Increases", "Decreases"),
            `Effect (log2FC)` = log2FoldChange,
            `SE (log2FC)` = lfcSE,
            `LRT Statistic` = stat,
            `p-value` = pvalue,
            `q-value` = qval,
            Signif = if_else(`q-value` < 0.05, "*", "")) %>%
  mutate_if(is.numeric, format, digits = 2)

# male samples
DEresults_male <-
  resultsmod2_malelinest %>%
  arrange(pvalue) %>%
  filter(row %in% DEresults_interaction$row) %>%
  left_join(annotation_miRNAs, by = c("row" = "miRNA"))

DEresults_male <-
  DEresults_male %>%
  filter(!Chr %in% c("chrX", "chrY"))

DEresults_male$qval <-
  qvalue(DEresults_male$pvalue)$qvalue

FinalResults_male <-
  DEresults_male %>%
  arrange(pvalue) %>%
  transmute(FC = 2^log2FoldChange,
            `Trajectory` = if_else(log2FoldChange > 0, "Increases", "Decreases"),
            miRNA = row,
            `Effect (log2FC)` = log2FoldChange,
            `SE (log2FC)` = lfcSE,
            `LRT Statistic` = stat,
            `p-value` = pvalue,
            `q-value` = qval,
            Signif = if_else(`q-value` < 0.05, "*", "")) %>%
  mutate_if(is.numeric, format, digits = 2) %>%
  filter(Signif == "*")


names(Table_S3)[-1] <-
  names(Table_S3)[-1] %>%
  paste0(., " [Interaction]")




# final merged table for reporting ---------------------------------------------
# ("addtraj": adding male and female trajectory information)

Table_S3_addtraj <-
FinalResults_male %>%
  left_join(FinalResults_female,
            by = "miRNA",
            suffix = c(" [Male]", " [Female]")) %>%
  left_join(Table_S3,
            .,
            by = "miRNA",
            suffix = c("[Interaction]", ""))
  
