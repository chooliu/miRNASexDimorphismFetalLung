# ==============================================================================
# 17_mod1.R
# perform RUVr && generate DESeq2 main effects model (no interactions)
# ==============================================================================




# RUVr calculations ------------------------------------------------------------

preruv_formula_model1 <-
  as.formula("~ AgeLin + AgeSq + Sex + SmokeBinary + run")
phenodata_final_model1 <-
  phenodata_final

vst_counts_prelim_model1 <-
  DESeqDataSetFromMatrix(
    countData = counts_final,
    design = preruv_formula_model1,
    colData = phenodata_final_model1
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions(fitType = "local") %>%
  varianceStabilizingTransformation(fitType = "local")

vst_residuals_model1 <-
  vst_counts_prelim_model1 %>%
  assay %>%
  apply(., 1,
        function(y) {
          tmp <-
            phenodata_final_model1 %>%
            bind_cols(y = y)
          lm(update(preruv_formula_model1, y ~ .),
             tmp) %>% 
            resid() }
)

check_ruv_kvals <-
  function(k) {
    ruv_tmp <-
      RUVr(x = assay(vst_counts_prelim_model1),
           k = k,
           residuals = vst_residuals_model1 %>% t,
           isLog = T)$W
    adonis(
      vst_residuals_model1 ~ ruv_tmp,
      method = "manhattan") %>%
      .$aov.tab %>%
      .$R2 %>%
      .[1]
  }

# choose k = 4 based on elbow plot
r2_by_k <- sapply(1:20, check_ruv_kvals)
ggplot(data = NULL, aes(x = 1:20, y = r2_by_k)) +
  geom_point() +
  geom_line() +
  xlab("# of RUVr Components (k)") +
  ylab("% of Residual Variance Explained\n(Multivariate PERMANOVA)") +
  theme_classic()

ruv_soln_model1 <-
  RUVr(x = assay(vst_counts_prelim_model1),
       k = 4,
       residuals = vst_residuals_model1 %>% t,
       isLog = T)

phenodata_final_model1 <-
  phenodata_final %>%
  bind_cols(
    ruv_soln_model1$W %>%
      as_tibble() %>%
      set_names(paste0("RUV", 1:4)))




# DESeq2 model (#1) ------------------------------------------------------------

base_formula_model1 <-
  update(preruv_formula_model1,  ~ . + RUV1 + RUV2 + RUV3 + RUV4)

deseqiterations <- 10000

deseq_obj_model1 <-
  DESeqDataSetFromMatrix(
    countData = 
      counts_final,
    design = base_formula_model1,
    colData = phenodata_final_model1
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions(fitType = "local") %>%
  nbinomWaldTest(maxit = deseqiterations)




# statistical testing ----------------------------------------------------------

# effect of nuisance variables, for reference
resultsmod1_run <-
  nbinomLRT(deseq_obj_model1,
            full = base_formula_model1,
            reduced = update(base_formula_model1, ~ . - run),
            maxit = deseqiterations)

resultsmod1_ruv <-
  nbinomLRT(deseq_obj_model1,
            full = base_formula_model1,
            reduced = preruv_formula_model1,
            maxit = deseqiterations)

# age
resultsmod1_age <-
  nbinomLRT(deseq_obj_model1,
            full = ~ Sex + SmokeBinary + run + RUV1 + RUV2 + RUV3 + 
              RUV4 + AgeSq + AgeLin,
            reduced = ~ Sex + SmokeBinary + run + RUV1 + RUV2 + RUV3 + RUV4 ,
            maxit = deseqiterations)

resultsmod1_agelin <-
  nbinomLRT(deseq_obj_model1,
            full =  ~ Sex + SmokeBinary + run + RUV1 + RUV2 + RUV3 + 
              RUV4 + AgeSq + AgeLin,
            reduced = ~ Sex + SmokeBinary + run + RUV1 + RUV2 + RUV3 + 
              RUV4 + AgeSq,
            maxit = deseqiterations)

resultsmod1_agesq <-
  nbinomLRT(deseq_obj_model1,
            full =  ~ Sex + SmokeBinary + run + RUV1 + RUV2 + RUV3 + 
              RUV4 + AgeLin + AgeSq,
            reduced = ~ Sex + SmokeBinary + run + RUV1 + RUV2 + RUV3 + 
              RUV4 + AgeLin,
            maxit = deseqiterations)

# sex, IUS
resultsmod1_sex <-
  nbinomLRT(deseq_obj_model1,
            full = ~ AgeLin + AgeSq + SmokeBinary + run + RUV1 + RUV2 + RUV3 + 
              RUV4 + Sex,
            reduced = ~ AgeLin + AgeSq + SmokeBinary + run +
              RUV1 + RUV2 + RUV3 + RUV4,
            maxit = deseqiterations)

resultsmod1_IUS <-
  nbinomLRT(deseq_obj_model1,
            full = ~ AgeLin + AgeSq + Sex + run + RUV1 + RUV2 + RUV3 + 
              RUV4 + SmokeBinary,
            reduced = ~ AgeLin + AgeSq + Sex + run +
              RUV1 + RUV2 + RUV3 + RUV4,
            maxit = deseqiterations)




# PCA visualizations (re-do now that duplicates are removed) -------------------

pca_soln <-
  ruv_soln_model1$normalizedCounts %>%
  t %>%
  prcomp(center = T, scale = F) # make sure this is subject x gene

pca_varexp <-
  (pca_soln$sdev^2 / sum(pca_soln$sdev^2) * 100) %>%
  head %>%
  format(digits = 2) %>%
  paste0(., "%")


phenodata_fmt_model1 <-
  phenodata_final %>%
  bind_cols(.,
            PC1 = pca_soln$x[ , 1],
            PC2 = pca_soln$x[ , 2]) %>%
  mutate(Sex = if_else(Sex, "Male", "Female"),
         IUS = case_when(SmokeBinary == 1 ~ "Yes",
                         SmokeBinary == 0 ~ "No",
                         ! SmokeBinary %in% c(1, 0) ~ "Unknown"),
         LibSize = log10(totcounts + 1),
         Run = factor(run, levels = unique(phenodata_final$run)[c(4, 1:3, 5:7)]),
         `Re-Run` = if_else(grepl('Rerun', run), "yes", "no"),
         RIN = as.numeric(RIN)) %>%
  left_join(phenodata_filtered %>% select(ALIASID, YearReceived),
            by = c("Recruitment_ID" = "ALIASID"))

phenodata_fmt_model1$totcounts <- assay(deseq_obj_model1) %>% colSums
phenodata_fmt_model1$miRNAs_Detected <- (assay(deseq_obj_model1) != 0) %>% colSums

levels(phenodata_fmt_model1$Run) <- 1:7

plot_pca_continuous <-
  function(variable_of_interest) {
    ggplot(data = phenodata_fmt_model1,
           aes_string("PC1", "PC2", color = variable_of_interest)) + 
      geom_point(alpha = 0.5, size = 1.25) +
      theme_classic() +
      xlab(paste0("PC1 (", pca_varexp[1], ")")) +
      ylab(paste0("PC2 (", pca_varexp[2], ")")) +
      scale_color_viridis(option = "B", direction = -1) +
      scale_fill_viridis(option = "B", direction = -1)
  }


plot_pca_categorical <-
  function(variable_of_interest) {
    ggplot(data = phenodata_fmt_model1,
           aes_string("PC1", "PC2",
                      shape = variable_of_interest,
                      color = variable_of_interest)) + 
      geom_point(alpha = 0.5, size = 1.25) +
      theme_classic() +
      xlab(paste0("PC1 (", pca_varexp[1], ")")) +
      ylab(paste0("PC2 (", pca_varexp[2], ")")) +
      scale_color_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      scale_shape_manual(values = c(0, 2, 3, 15, 16, 17, 18))
  }


plot_grid(
  plot_pca_continuous("Age"),
  plot_pca_categorical("IUS"),
  plot_pca_categorical("Sex") + scale_color_brewer(palette = "Set1"),
  plot_pca_categorical("Superpopulation") + scale_color_brewer(palette = "Set2"),
  plot_pca_continuous("genoPC1") + scale_color_viridis(option = "C"),
  plot_pca_continuous("genoPC2") + scale_color_viridis(option = "C"),
  ncol = 2
)

plot_grid(
  plot_pca_continuous("RIN") + scale_color_viridis(option = "A"),
  plot_pca_continuous("LibSize") + scale_color_viridis(option = "A", direction = -1),
  plot_pca_continuous("miRNAs_Detected") + ylab("# miRNAs"),
  nrow = 2
)





