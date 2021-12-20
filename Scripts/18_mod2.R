# ==============================================================================
# 18_mod2.R
# repeat RUV & DESeq2 in script 17, but with Age*Sex interaction
# ==============================================================================




# repeat RUVr calculations -----------------------------------------------------

preruv_formula_model2 <-
  as.formula("~ AgeLin + AgeSq + Sex + AgeLin*Sex + SmokeBinary + run")
phenodata_final_model2 <- phenodata_final

vst_counts_prelim_model2 <-
  DESeqDataSetFromMatrix(
    countData = counts_final,
    design = preruv_formula_model2,
    colData = phenodata_final_model2
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions(fitType = "local") %>%
  varianceStabilizingTransformation(fitType = "local")

vst_residuals_model2 <-
  vst_counts_prelim_model2 %>%
  assay %>%
  apply(., 1,
        function(y) {
          tmp <-
            phenodata_final_model2 %>%
            bind_cols(y = y)
          lm(update(preruv_formula_model2, y ~ .),
             tmp) %>% 
            resid() }
  )

check_ruv_kvals <-
  function(k) {
    ruv_tmp <-
      RUVr(x = assay(vst_counts_prelim_model2),
           k = k,
           residuals = vst_residuals_model2 %>% t,
           isLog = T)$W
    adonis(
      vst_residuals_model2 ~ ruv_tmp,
      method = "manhattan") %>%
      .$aov.tab %>%
      .$R2 %>%
      .[1]
  }


# choose k = 4
r2_by_k <- sapply(1:20, check_ruv_kvals)
ggplot(data = NULL, aes(x = 1:20, y = r2_by_k)) +
  geom_point() +
  geom_line() +
  xlab("# of RUVr Components (k)") +
  ylab("% of Residual Variance Explained\n(Multivariate PERMANOVA)") +
  theme_classic()

ruv_soln_model2 <-
  RUVr(x = assay(vst_counts_prelim_model2),
       k = 4,
       residuals = vst_residuals_model2 %>% t,
       isLog = T)

phenodata_final_model2 <-
  phenodata_final %>%
  bind_cols(
    ruv_soln_model2$W %>%
      as_tibble() %>%
      set_names(paste0("RUV", 1:4)))




# DESeq2 model (#2) ------------------------------------------------------------

base_formula_model2 <-
  update(preruv_formula_model2,  ~ . + RUV1 + RUV2 + RUV3 + RUV4)

deseq_obj_model2 <-
  DESeqDataSetFromMatrix(
    countData = 
      counts_final,
    design = base_formula_model2,
    colData = phenodata_final_model2
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions(fitType = "local") %>%
  nbinomWaldTest(maxit = deseqiterations)


plotDispEsts(deseq_obj_model2)




# statistical testing ----------------------------------------------------------
# focused on sex, interact, sex-specific trajectories (Wald)

# interaction
resultsmod2_firstorderinteract <-
  nbinomLRT(deseq_obj_model2,
            full = base_formula_model2,
            reduced = update(base_formula_model2, ~ . - AgeLin:Sex),
            maxit = deseqiterations)

# male age trajectory
resultsmod2_malelinest <-
  results(deseq_obj_model2,
          contrast = as.numeric(colnames(coef(deseq_obj_model2)) %in%
                                  c("AgeLin", "AgeLin.SexTRUE")),
          tidy = T)

# female age trajectory
resultsmod2_femalelinest <-
  results(deseq_obj_model2,
          contrast = as.numeric(colnames(coef(deseq_obj_model2)) %in%
                                  c("AgeLin")),
          tidy = T)

# sex primary effect, but with model 2 
resultsmod2_sex <-
  results(deseq_obj_model2,
          contrast = as.numeric(colnames(coef(deseq_obj_model2)) %in%
                                  c("SexTRUE")),
          tidy = T)
