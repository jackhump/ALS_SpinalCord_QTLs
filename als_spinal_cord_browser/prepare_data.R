# prepare shared data for visualisation
library(tidyverse)

inFolder <- "../NYGC_ALS/ALS_SC/sharing"
#gene_meta <- read_tsv()
tissues <- c("Cervical_Spinal_Cord", "Lumbar_Spinal_Cord", "Thoracic_Spinal_Cord")

tpm_all <- map(tissues, ~{
  read_tsv(paste0(inFolder, "/", .x, "/", .x, "_gene_tpm.tsv.gz" )) }) %>%
  reduce(inner_join, by = c("ensembl_id", "gene_name"))

metadata <- map_df(tissues, ~{
  read_tsv(paste0(inFolder, "/", .x, "/", .x, "_metadata.tsv.gz" ))
}) %>%
  filter(disease %in% c("Control", "ALS"))

counts_all <- map(tissues, ~{
  read_tsv(paste0(inFolder, "/", .x, "/", .x, "_gene_counts.tsv.gz" )) }) %>%
  reduce(inner_join, by = c("ensembl_id", "gene_name"))

# test matching - only males should have UTY expression
filter(tpm_all, gene_name == "UTY") %>% 
  pivot_longer(names_to = "rna_id", values_to = "TPM", cols = !c(ensembl_id, gene_name)) %>%
  left_join(metadata, by = c("rna_id") ) %>%
  ggplot(aes(x = sex, y = TPM)) + geom_boxplot() + geom_point()

filter(counts_all, gene_name == "UTY") %>% 
  pivot_longer(names_to = "rna_id", values_to = "count", cols = !c(ensembl_id, gene_name)) %>%
  left_join(metadata, by = c("rna_id") ) %>%
  ggplot(aes(x = sex, y = count)) + geom_boxplot() + geom_point()

#tpm_all <- tpm_all[ ,metadata$rna_id]
#counts_all <- counts_all[, metadata$rna_id]

keep_exp <- rowSums(tpm_all > 0 ) > 0.5 * ncol(tpm_all)


prep_matrix <- function(df, keep){
  df <- df[ keep, c("gene_name", metadata$rna_id)]
  df <- df[ !duplicated(df$gene_name) & !is.na(df$gene_name),]
  df <- as.data.frame(df)
  row.names(df) <- df$gene_name
  df$ensembl_id <- NULL
  df$gene_name <- NULL
  return(df)
}

tpm_df <- prep_matrix(tpm_all, keep = keep_exp)
counts_df <- prep_matrix(counts_all, keep_exp)

# check again
as.data.frame(tpm_df[ "UTY",]) %>% rownames_to_column("gene") %>%
  pivot_longer(names_to = "rna_id", values_to = "count", cols = !c(gene)) %>%
  left_join(metadata, by = "rna_id") %>%
  ggplot(aes(x = sex, y = count)) + geom_boxplot() + geom_point()
as.data.frame(counts_df[ "UTY",]) %>% rownames_to_column("gene") %>%
  pivot_longer(names_to = "rna_id", values_to = "count", cols = !c(gene)) %>%
  left_join(metadata, by = "rna_id") %>%
  ggplot(aes(x = sex, y = count)) + geom_boxplot() + geom_point()


## METADATA

metadata <- metadata %>%
  mutate(disease_c9 = case_when(
    disease == "Control" ~ "Control",
    disease == "ALS" & mutations == "C9orf72" ~ "ALS-C9",
    disease == "ALS" & mutations == "None" ~ "ALS",
    TRUE ~ "ALS-Other"
  ))


# metadata <- metadata %>%
#   #mutate(sample = external_sample_id) %>%
#   select(-age_at_death, -external_sample_id, -external_subject_id, -site_specimen_collected, -differential_expression, -qtl_mapping)

metadata$tissue <- gsub("_Spinal_Cord", "", metadata$tissue)

metadata$disease = factor(metadata$disease, levels = c("Control", "ALS"))
metadata$`C9orf72 status` <- factor(metadata$disease_c9, levels = c("Control", "ALS", "ALS-C9", "ALS-Other"))

metadata$RIN <- metadata$rin
metadata$`disease duration` <- metadata$disease_duration
metadata$site_of_motor_onset <- gsub(" ", "\n",metadata$site_of_motor_onset )

## add in de and duration results

load("../../NYGC_ALS/data/differential_expression/ALS/ALS_Spinal_Cord_de_res_Mar2022.RData")
de_res <- de_res_final
de_res$TSC <- NULL
names(de_res) <- c("Cervical", "Lumbar")

dur_res <- load("../../NYGC_ALS/data/differential_expression/ALS/Model_8_limma_res_Mar2022.RData")
dur_res <- de_res_final
dur_res$TSC <- NULL
names(dur_res) <- c("Cervical", "Lumbar")



save(metadata, tpm_df, de_res, dur_res, file = here::here("tpm_metadata.RData"))



stop()

## TMM, VOOM and RBE
library(limma)
library(edgeR)
gExpr <- DGEList(counts_df)
# Perform TMM normalization
gExpr <- calcNormFactors(gExpr,method = "TMM")
# Specify variables to be included in the voom() estimates of
# uncertainty.
# Recommend including variables with a small number of categories
# that explain a substantial amount of variation
design <- model.matrix( ~ 1, metadata)

counts_voom <-voom(gExpr, design = design)$E

# removeBatchEffects
covariate_df <- 
   metadata %>%
   dplyr::select(RIN, library_prep) %>%
   mutate(library_prep = as.numeric(as.factor(library_prep)) )
# 
row.names(covariate_df) <- metadata$sample
# 

counts_voom_rbe <- removeBatchEffect(counts_voom, batch = metadata$site_id, covariates = covariate_df )

metadata$sample <- metadata$rna_id

counts_voom  <- as.data.frame(counts_voom)

save(metadata, counts_voom, file = here::here("../als_spinal_cord_browser/voom_metadata.RData"))

# CHECK XIST AND SEX - WEIRD!
