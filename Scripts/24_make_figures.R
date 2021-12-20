# ==============================================================================
# 24_make_figures.R
# generate main text figures (Fig1, 2, 3)
# ==============================================================================




# Fig 1 ------------------------------------------------------------------------
# check distribution of age, by sex

wilcox.test(Age ~ Sex, phenodata_final_model1)

ggplot(phenodata_final_model1,
       aes(Age, Sex, fill = Sex)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 1, height = 0, seed = 1),
    point_shape = '|', point_size = 3, point_alpha = 1,
    alpha = 0, scale = 1.5, lty = 0) +
  geom_density_ridges(aes(height = ..density..),
                      alpha = 0.4, scale = 1.5,
                      stat = "density", trim = T) +
  xlab("Estimated Gestational Age (days)") +
  scale_y_discrete(labels = c("Male", "Female"),
                   expand = expansion(mult = c(0.01, .7))) +
  theme_ridges(font_size = 12, grid = F)  +
  scale_fill_manual(values = c("maroon", "black")) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(50, 120), breaks = c(54, 70, 100, 117)) +
  coord_cartesian(clip = "off")

ggsave("./FiguresTables/Fig1.png", width = 4, height = 2.5)





# Fig 2 ------------------------------------------------------------------------
# 2a: "forest plot" showing top 10 miRNAs with 95% CIs &
# 2bc: examples of miRNAs diff exp by sex (main effect), adj for covariates

# 2a
feature_miRNA_covar1 <- 
  apply(feature_miRNA_vst, 1, 
        function(y) {
          lm(y ~ 
               run + RUV1 + RUV2 + RUV3 + RUV4, phenodata_final_model1) %>% resid
        })

resultssex_top10 <-
  resultsmod1_sex %>%
  results(tidy = T) %>% 
  filter(row %in% FinalResults_Sex$miRNA[1:10]) %>%
  mutate(mid = log2FoldChange,
         low = mid - 1.96*lfcSE,
         high = mid + 1.96*lfcSE,
         dir = mid > 0) %>%
  arrange(mid) %>%
  left_join(., FinalResults_Sex, by = c(row = "miRNA")) %>%
  mutate(miRNA = row %>% as.factor %>% fct_inorder())

fig2a_bot <-
  ggplot(data = resultssex_top10) +
  geom_errorbar(aes(x = miRNA, ymin = low, ymax = high,
                    color = dir), width = 0.4, size = 0.5) +
  geom_point(aes(x = miRNA, y = mid, color = dir), shape = 15) +
  geom_hline(yintercept = 0, lty = 3) +
  coord_flip() +
  theme_few() +
  xlab("miRNA") +
  scale_y_continuous(expression(log[2]~"(Fold Change)"),
                     limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "maroon"))

fig2a_top <-
  ggplot(data = NULL) +
  annotate(geom = "text", label = "higher in female",
           x = -2, y  = 1, hjust = 0, color = "black", size = 3) +
  annotate(geom = "text", label = "higher in male",
           x = 5, y = 1, hjust = 0, color = "maroon", size = 3) +
  annotate("segment", x = 3.5, xend = -0.5, y = 0.25, yend = 0.25, size = 0.25,
           arrow = arrow(type = "closed", angle = 20, length = unit(0.3, units = "cm")),
           color = "black") +
  annotate("segment", x = 5, xend = 9, y = 0.25, yend = 0.25, size = 0.25,
           arrow = arrow(type = "closed", angle = 20, length = unit(0.3, units = "cm")),
           color = "maroon") +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_void()

fig2a_right <-
  ggplot(data = resultssex_top10 %>%
           mutate(FC = 2^mid %>% formatC(., digits = 2, format = "f", flag = "#"),
                  pformat = gsub("e-", "%*%10^{-", `q-value`) %>% paste0(., "}") %>%
                    as.expression %>% as.character),
         aes(y = miRNA, x = 0)) +
  geom_text(aes(label = FC, x = 0), hjust = 0, size = 2.75) +
  geom_text(aes(label = pformat, x = 0.5), hjust = 0, parse  = T, size = 2.75) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_void()

# 2bc
plot_mirna_by_sex <-
  function(name_of_mirna, residuals = T) {
    
    if (residuals) {
      tmp_plot_df <-
        tibble(Sex = if_else(phenodata_final_model2$Sex, "Male", "Female"),
               Age = phenodata_final_model2$Age,
               y = feature_miRNA_covar1 %>% .[, name_of_mirna])
      ggplot(tmp_plot_df, aes(Sex, y, color = Sex)) + 
        geom_quasirandom(alpha = 0.5, bandwidth = 0.2) +
        theme_few() +
        xlab("Sex") + ylab(name_of_mirna) +
        stat_summary(geom = "errorbar", fun.min = mean, color = "black") +
        scale_color_manual(values = palette_color_sex) +
        scale_shape_manual(values = palette_shape_sex) +
        theme(legend.position = "none")
    } else {
      tmp_plot_df <-
        tibble(Sex = if_else(phenodata_final_nodupes$Sex, "Male", "Female"),
               y = feature_miRNA_vst %>% assay %>% .[name_of_mirna, ])
      ggplot(tmp_plot_df, aes(Sex, y, color = Sex)) + 
        geom_quasirandom(alpha = 0.5, bandwidth = 0.25) +
        theme_few() +
        stat_summary(geom = "errorbar", fun = mean) +
        xlab("Sex") + ylab(name_of_mirna) +
        scale_color_manual(values = palette_color_sex) +
        scale_shape_manual(values = palette_shape_sex) +
        theme(legend.position = "none")
    }
  }

fig2bc <-
  plot_grid(
    plot_mirna_by_sex("hsa-miR-27b-3p", T),
    plot_mirna_by_sex("hsa-miR-31-5p", T),
    labels = c("(B)", "(C)"), label_size = 10)

# putting it all together
fig2 <-
  plot_grid(
    plot_grid(fig2a_top, fig2a_bot,
              rel_heights = c(1, 9),
              nrow = 2),
    plot_grid(ggplot(data = NULL) + theme_void(),
              fig2a_right, ggplot(data = NULL) + theme_void(),
              rel_heights = c(1.5, 9, 2.25), nrow = 3),
    rel_widths = c(9, 3))

plot_grid(fig2,
          fig2bc,
          nrow = 2,
          rel_heights = c(4, 3),
          labels = c("(A)", ""), label_size = 10)
ggsave(filename = "Rev_Fig2.pdf", width = 5, height = 6)




# Fig 3 ------------------------------------------------------------------------
# examples of sex-specific age-trajectories (interaction; adj for covariates)

feature_miRNA_covar2 <-
  apply(feature_miRNA_vst, 1,
        function(y) {
          lm(y ~ AgeLin + AgeSq + SmokeBinary + run +
               RUV1 + RUV2 + RUV3 + RUV4, phenodata_final_model2) %>% resid
        })

plot_mirna_by_agesex <-
  function(name_of_mirna, residuals = T) {
    
    if (residuals) {
      tmp_plot_df <-
        tibble(Sex = if_else(phenodata_final_model2$Sex, "Male", "Female"),
               Age = phenodata_final_model2$Age,
               y = feature_miRNA_covar2 %>% .[, name_of_mirna])
      ggplot(tmp_plot_df, aes(Age, y, color = Sex, shape = Sex, lty = Sex)) + 
        geom_quasirandom(alpha = 0.2, bandwidth = 0.2) +
        theme_few() +
        xlab("Predicted Age (days)") + scale_y_continuous(name_of_mirna, breaks = pretty_breaks(n = 5)) +
        geom_smooth(method = "lm", se = F) +
        scale_color_manual(values = c("black", "maroon")) +
        scale_shape_manual(values = c(15, 19)) +
        scale_linetype_manual(values = c(1, 2))
    } else {
      tmp_plot_df <-
        tibble(Sex = if_else(phenodata_final_nodupes$Sex, "Male", "Female"),
               y = feature_miRNA_vst %>% assay %>% .[name_of_mirna, ])
      ggplot(tmp_plot_df, aes(Age, y, color = Sex, shape = Sex, lty = Sex)) + 
        geom_quasirandom(alpha = 0.2, bandwidth = 0.25) +
        theme_few() +
        geom_smooth(method = "lm", se = F) +
        xlab("Sex") + scale_y_continuous(name_of_mirna, breaks = pretty_breaks(n = 5)) +
        scale_color_manual(values = palette_color_sex) +
        scale_shape_manual(values = palette_shape_sex) +
        scale_linetype_manual(values = c(1, 2))
    }
  }

fig3 <-
  plot_grid(
    plot_mirna_by_agesex("hsa-let-7b-5p", T) + theme(legend.position = "none"),
    plot_mirna_by_agesex("hsa-miR-27b-3p", T) + theme(legend.key.width = unit(1, "cm")),
    rel_widths = c(4, 6), 
    labels = c("(A)", "(B)"), label_size = 10)

ggsave(fig3,filename =  "./FiguresTables/Fig3.png",
       height = 3, width = 7)
system('open  "./FiguresTables/Fig3.png"')
