# ==============================================================================
# 23_interact_model_sensanalysis.R 
# script 21, but check sex main effects in interaction model (script 18)
# ==============================================================================



# repeat calculations in script 21  --------------------------------------------
# but based on "alt"-ernative model (model 2 / interactive model)

DEresults_sexalt <-
  resultsmod2_sex %>%
  arrange(pvalue) %>%
  filter(!is.na(padj)) %>%
  left_join(annotation_miRNAs %>% filter(!duplicated(miRNA)),
            by = c("row" = "miRNA"))

DEresults_sexalt <-
  DEresults_sexalt %>%
  filter(!Chr %in% c("chrX", "chrY"))

DEresults_sexalt$qval <-
  qvalue(DEresults_sexalt$pvalue)$qvalue

FinalResults_sexalt <-
  DEresults_sexalt %>%
  arrange(pvalue) %>%
  rowwise() %>%
  transmute(FC = 2^log2FoldChange,
            `Higher Levels In` = if_else(log2FoldChange > 0, "Males", "Females"),
            miRNA = row,
            `Effect (log2FC)` = log2FoldChange,
            `SE (log2FC)` = lfcSE,
            `LRT Statistic` = stat,
            `p-value` = pvalue,
            `q-value` = qval,
            Signif = if_else(`q-value` < 0.05, "*", "")) %>%
  mutate_if(is.numeric, format, digits = 2) 

Table_S1alt <-
  FinalResults_sexalt %>%
  filter(Signif == "*")




# compare miRNAs identified as DE, btwn models ---------------------------------
# sensitivity analysis suggests very similar results / high concordance
# below figures presented in article supplementary material

# identity of signif miRNAs
venn(list(`Main Effects` = Table_S1$miRNA,
          `Interactive Only` = Table_S1alt$miRNA)) %>%
  plot

setdiff(Table_S1$miRNA, Table_S1alt$miRNA)
setdiff(Table_S1alt$miRNA, Table_S1$miRNA)


# similarity in log2FCs
left_join(DEresults_sex, DEresults_sexalt,
          by = "row") %>%
  select(contains("log2Fold")) %>%
  cor(use = "complete")

left_join(DEresults_sex, DEresults_sexalt,
          by = "row") %>% 
  ggplot(data = ., aes(log2FoldChange.x, log2FoldChange.y)) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_few() +
  xlab("log2(Fold Change)\nMain-Effect Only Model") +
  ylab("log2(Fold Change)\nInteractive-Effects Model") +
  annotate(x = -0.8, y = 0.6, geom = "text", label = "Pearson r = 0.999")


# similarity in p-values
left_join(DEresults_sex, DEresults_sexalt,
          by = "row") %>%
  select(contains("pval")) %>%
  mutate_all(log10) %>%
  cor(use = "complete")

left_join(DEresults_sex, DEresults_sexalt,
          by = "row") %>% 
  ggplot(data = ., aes(pvalue.x, pvalue.y)) +
  geom_point(alpha = 0.25) +
  theme_few() +
  scale_x_log10("p-value (log10-Scale)\nMain-Effect Only Model") +
  scale_y_log10("p-value (log10-Scale)\nInteractive-Effects Model") +
  annotate(x = 5e-9, y = 1, geom = "text",
           label = "Pearson r = 0.999")




