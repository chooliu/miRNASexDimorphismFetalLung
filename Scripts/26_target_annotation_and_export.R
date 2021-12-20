# ==============================================================================
# 26_target_annotation_and_export.R
# annotate results tables with putative mRNA targets
# ==============================================================================




# load libraries for miRNAtap --------------------------------------------------
# not incl. script 01 due to some function names / scoping issues
# (namely, Dbi select versus tidyverse select)

library(biomaRt)
library(miRNAtap)


# biomart for entrez --> symbol name -------------------------------------------
# (particular ensembl release frozen to match some other analyses underway)

biomart_obj <-
  useMart(biomart = "ENSEMBL_MART_ENSEMBL",
          dataset = "hsapiens_gene_ensembl",
          host = "sep2019.archive.ensembl.org") %>%
  getBM(attributes =
          c("refseq_mrna", "external_gene_name",
            "ensembl_gene_id", "ensembl_transcript_id", "entrezgene_id"),
        mart = .)



# helper fxn to query & fmt miRNAtap db ----------------------------------------

get_predicted_symbols <- 
  function(mirname) {
    
  predicted_out <-
    getPredictedTargets(mirname, species = 'hsa')
    
    predicted_out %>%
      row.names() %>%
      .[1:min(length(.), 10)] %>%
      tibble(entrezgene_id = as.integer(.)) %>%
      left_join(biomart_obj) %>%
      group_by(external_gene_name, entrezgene_id) %>%
      group_keys() %>%
      .$external_gene_name %>%
      paste0(collapse = "; ")
  }



# run above fxn to add predicted targets to results tables ---------------------

Table_S1 <-
  Table_S1 %>%
  rowwise() %>%
  mutate(Predicted = get_predicted_symbols(miRNA))


Table_S3_addtraj <-
  Table_S3_addtraj %>%
  rowwise() %>%
  mutate(Predicted = get_predicted_symbols(miRNA))



# export final results tables --------------------------------------------------

list(Table_S1 = Table_S1) %>%
  write_xlsx("./FiguresTables/TableS1.xlsx", format_headers = F)

list(Table_S3 = Table_S3_addtraj) %>%
  write_xlsx("./FiguresTables/TableS3.xlsx", format_headers = F)

