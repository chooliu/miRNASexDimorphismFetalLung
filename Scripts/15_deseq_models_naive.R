# ==============================================================================
# 15_deseq_models_naive.R
# "naive" model: no RUVSeq, has duplicated samples (tech reps/re-sequence)
# for initial visualization esp. assessment of batch effects
# ==============================================================================




# generate "naive" DESeq2 model ------------------------------------------------

preruv_formula <-
  as.formula("~ AgeLin + AgeSq + Sex + SmokeBinary + run")

deseq_obj_naive <-
  DESeqDataSetFromMatrix(
    countData = 
      counts_naive %>%
      select(-run_sample, -Recruitment_ID) %>%
      t(),
    design = preruv_formula,
    colData = phenodata_naive
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions(fitType = "local")

plotDispEsts(deseq_obj_naive)

vst_counts_naive <-
  deseq_obj_naive %>%
  varianceStabilizingTransformation(fitType = "local")




# PCA solution and visualization -----------------------------------------------

# run PCA
pca_soln_naive <-
  assay(vst_counts_naive) %>%
  t %>%
  prcomp(center = T, scale = F) # make sure this is subject x gene

pca_varexp_naive <-
  (pca_soln_naive$sdev^2 / sum(pca_soln_naive$sdev^2) * 100) %>%
  head %>%
  format(digits = 2) %>%
  paste0(., "%")

# better formatted metadata for plotting
# note: run ordering is due to alphabetical order != date run
phenodata_naive_fmt <-
  phenodata_naive %>%
  bind_cols(.,
            PC1 = pca_soln_naive$x[ , 1],
            PC2 = pca_soln_naive$x[ , 2]) %>%
  mutate(Sex = if_else(Sex, "Male", "Female"),
         IUS = case_when(SmokeBinary == 1 ~ "Yes",
                         SmokeBinary == 0 ~ "No",
                         ! SmokeBinary %in% c(1, 0) ~ "Unknown"),
         LogLibSize = log10(totcounts + 1),
         Run = factor(run,
                      levels = unique(phenodata_naive$run)[c(4, 1:3, 5:7)]),
         `Re-Run` = if_else(grepl('Rerun', run), "yes", "no"),
         RIN = as.numeric(RIN)) %>%
  left_join(phenodata_filtered %>% select(ALIASID, YearReceived),
            by = c("Recruitment_ID" = "ALIASID"))


# add total counts for visualization
phenodata_naive_fmt$totcounts <- assay(deseq_obj_naive) %>% colSums
phenodata_naive_fmt$miRNAs_Detected <- (assay(deseq_obj_naive) != 0) %>% colSums

# tidy run names
levels(phenodata_naive_fmt$Run) <- 1:7

# functions to plot variables in "phenodata_naive_fmt" quickly
plot_naive_pca_continuous <-
  function(variable_of_interest) {
    ggplot(data = phenodata_naive_fmt,
           aes_string("PC1", "PC2", color = variable_of_interest)) +
      geom_point(alpha = 0.5, size = 1.25) +
      theme_classic() +
      xlab(paste0("PC1 (", pca_varexp_naive[1], ")")) +
      ylab(paste0("PC2 (", pca_varexp_naive[2], ")")) +
      scale_color_viridis(option = "B", direction = -1) +
      scale_fill_viridis(option = "B", direction = -1) +
      theme(legend.text = element_text(size = 7),
            legend.title =  element_text(size = 10))
  }

plot_naive_pca_categorical <-
  function(variable_of_interest) {
      ggplot(data = phenodata_naive_fmt, aes_string("PC1", "PC2",
                                  shape = variable_of_interest,
                                  color = variable_of_interest)) +
      geom_point(alpha = 0.5, size = 1.25) +
      theme_classic() +
      xlab(paste0("PC1 (", pca_varexp_naive[1], ")")) +
      ylab(paste0("PC2 (", pca_varexp_naive[2], ")")) +
      scale_color_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      scale_shape_manual(values = c(0, 2, 3, 15, 16, 17, 18)) +
      theme(legend.text = element_text(size = 7),
            legend.title =  element_text(size = 10))
  }

# examine PCA, colored by known biological features
plot_grid(
  plot_naive_pca_continuous("Age"),
  plot_naive_pca_categorical("IUS"),
  plot_naive_pca_categorical("Sex") + scale_color_brewer(palette = "Set1"),
  plot_naive_pca_categorical("Superpopulation") + scale_color_brewer(palette = "Set2"),
  plot_naive_pca_continuous("genoPC1") + scale_color_viridis(option = "C"),
  plot_naive_pca_continuous("genoPC2") + scale_color_viridis(option = "C"),
  ncol = 2
)

# examine known technical features
plot_grid(
  plot_naive_pca_categorical("`Re-Run`"),
  plot_naive_pca_continuous("RIN") +
    scale_color_viridis(option = "A"),
  plot_naive_pca_continuous("LogLibSize") +
    scale_color_viridis(option = "A", direction = -1),
  plot_naive_pca_continuous("miRNAs_Detected"),
  nrow = 2
)




# check if replicates cluster near one another
# ------------------------------------------------------------------------------
set.seed(1111)
df_annotate_pca_replicates <-
  phenodata_naive_fmt %>%
  mutate(Replicates =
           if_else(Recruitment_ID %in% names(min_corrs_in_dupes), # [*]
                   Recruitment_ID, "") %>%
             as.factor()
  ) %>%
  filter(Replicates != "")
levels(df_annotate_pca_replicates$Replicates) <-
  LETTERS[1:length(levels(df_annotate_pca_replicates$Replicates))]

# [*] alt, replace with low or high correlation IDs only
# c(sample((min_corrs_in_dupes %>% .[`<`(., 0.95)] %>% names), 4),
#   sample((min_corrs_in_dupes %>% .[`>`(., 0.99)] %>% names), 4)),

phenodata_naive_fmt %>%
  ggplot(data = ., aes(PC1, PC2)) +
  geom_point(color = "grey", alpha = 0.5) +
  geom_text(data = df_annotate_pca_replicates,
            aes(PC1, PC2, label = Replicates, fill = NULL), size = 4) +
  geom_convexhull(data = df_annotate_pca_replicates,
                  aes(PC1, PC2, group = Replicates),
                  fill = NA, color = "black") +
  theme_classic() +
  xlab(paste0("PC1 (", pca_varexp_naive[1], ")")) +
  ylab(paste0("PC2 (", pca_varexp_naive[2], ")"))




# check associations between metadata features with PCA location, libsize
# ------------------------------------------------------------------------------
ggplot(phenodata_naive_fmt,
       aes(x = RIN, totcounts)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_smooth() +
  scale_y_log10("Log10(Total Counts)")
ggplot(phenodata_naive_fmt,
       aes(x = RIN, totcounts)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_smooth() +
  scale_y_log10("Log10(Total Counts)")

ggplot(phenodata_naive_fmt,
       aes(`Re-Run`, LogLibSize)) +
  geom_point()
ggplot(phenodata_naive_fmt,
       aes(run, LogLibSize)) +
  geom_point()




# DE by run --------------------------------------------------------------------
# check number of significantly DE miRNAs by run (299 of 899 have FDR < 0.05)
# note: other iterations did not adjust for run, but included explicitly for EDA
# (alternatives were to rely on the sample block-randomization by run /
# attempt to resolve run effects implicity via RUVSeq, but high # of run-assoc)

deseq_test_results_run <-
  nbinomLRT(deseq_obj_naive,
            full = ~ AgeLin + AgeSq + Sex + SmokeBinary + run,
            reduced = ~ AgeLin + AgeSq + Sex + SmokeBinary)

naive_results_by_run <-
  results(deseq_test_results_run, tidy = T) %>%
  filter(padj < 0.05) %>%
  arrange(padj)

nrow(naive_results_by_run)

plot_mirna_by_run <-
  function(name_of_mirna) {

    tmp_plot_df <-
      tibble(x = as.factor(phenodata_naive$run) %>%
               factor(., levels = levels(.)[c(4, 1:3, 5:7)]),
             y = vst_counts_naive %>% assay %>% .[name_of_mirna, ])

    ggplot(tmp_plot_df,
           aes(x, y, color = x)) +
      geom_quasirandom() +
      xlab("Run") + ylab(paste0(name_of_mirna, " (VST)")) +
      stat_summary(geom = "errorbar", fun.min = median, color = "black") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1),
            legend.position = "none")
  }

plot_grid(
  plot_mirna_by_run(naive_results_by_run$row[1]) +
  scale_x_discrete(labels = NULL),
  plot_mirna_by_run(naive_results_by_run$row[4]) +
  scale_x_discrete(labels = NULL),
  plot_mirna_by_run(naive_results_by_run$row[11]),
  nrow = 3,
  rel_heights = c(0.3, 0.3, 0.6))

