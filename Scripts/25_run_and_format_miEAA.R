# ==============================================================================
# 25_run_and_format_miEAA.R
# format results for export --> GSEA analysis (miEAA/miRWalk)
# ==============================================================================




# export list of miRNA ordered by p-value --------------------------------------
# text lists --> miEAA web portal; accessed Nov2021 --> export as spreadsheet
#  (settings: GSEA, miRWalk pathways, FDR < 0.10)

dir.create("miEAA_Input")
DEresults_sex$row %>%
  write_lines(file = "miEAA_Input/all_GSEA.txt")
DEresults_interaction$row %>%
  write_lines(file = "miEAA_Input/interact_GSEA.txt")




# take output of miEAA and select enriched pathways only  ----------------------
# (for ease of interpretation)

read_xlsx("./Data/miEAA/miEAA-Sex.xlsx", skip = 1) %>%
  filter(Enrichment == "enriched") %>%
  list(TableS2 = .) %>%
  write_xlsx("./FiguresTables/TableS2.xlsx", format_headers = F)

read_xlsx("./Data/miEAA/miEAA-Interaction.xlsx", skip = 1) %>%
  filter(Enrichment == "enriched") %>%
  list(TableS4 = .) %>%
  write_xlsx("./FiguresTables/TableS4.xlsx", format_headers = F)
