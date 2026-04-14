
# ============================================================
# PROJECT 1: DNA VARIANT ANALYSIS — COMPLETE PIPELINE
# Author: Stephen Opoku Asante
# Institution: University of Cape Coast, Ghana
# Date: April 2026
# ============================================================

rm(list = ls())
library(tidyverse)
library(vcfR)

project_folder <- "C:/Users/Stephen Opoku Asante/Desktop/Project_folder"
dir.create(project_folder, showWarnings = FALSE)
dir.create(paste0(project_folder, "/plots"), showWarnings = FALSE)
setwd(project_folder)

# SECTION 1: SIMULATE DATA
set.seed(42)
n_variants <- 500
n_samples  <- 20

samples <- c(paste0("PATIENT_", 1:10),
             paste0("CONTROL_", 1:10))
chromosomes <- paste0("chr", c(1:22))

vcf_data <- tibble(
  CHROM  = sample(chromosomes, n_variants, replace = TRUE),
  POS    = sample(1000000:250000000, n_variants, replace = FALSE),
  ID     = paste0("rs", sample(100000:999999, n_variants)),
  REF    = sample(c("A","T","C","G"), n_variants, replace = TRUE),
  ALT    = sample(c("A","T","C","G"), n_variants, replace = TRUE),
  QUAL   = round(runif(n_variants, min=20, max=60), 1),
  FILTER = sample(c("PASS","PASS","PASS","LowQual"),
                  n_variants, replace = TRUE)
) %>%
  filter(REF != ALT) %>%
  arrange(CHROM, POS)

# SECTION 2: QUALITY FILTERING
vcf_pass <- vcf_data %>% filter(FILTER == "PASS")

vcf_pass <- vcf_pass %>%
  mutate(
    SAMPLE = sample(samples, nrow(vcf_pass), replace = TRUE),
    GROUP  = ifelse(str_detect(SAMPLE, "PATIENT"),
                    "Patient", "Control"),
    DP     = round(rnorm(nrow(vcf_pass), mean = 30, sd = 6)),
    AF     = round(runif(nrow(vcf_pass), min = 0.25, max = 0.85), 3)
  ) %>%
  mutate(DP = ifelse(GROUP == "Patient", round(DP * 0.88), DP))

# SECTION 3: STATISTICAL TEST
qual_test <- t.test(QUAL ~ GROUP, data = vcf_pass)

# SECTION 4: ANNOTATION
disease_genes <- c("BRCA1","BRCA2","TP53","CFTR","LDLR",
                   "APOE","MTHFR","G6PD","HBB","HBA1","PCSK9")
background_genes <- c("ACTB","GAPDH","MYH9","TTN","MUC16",
                      "OBSCN","SYNE1","RYR2","DNAH5","HMCN1",
                      "FLG","ABCA13","LAMA1","CUBN","XIRP2",
                      "ADGRV1","USH2A","GPR98","AGRN","NEB")
consequence_types <- c("missense_variant","synonymous_variant",
                       "stop_gained","splice_region_variant",
                       "intron_variant")
clinical_significance <- c("Pathogenic","Likely_pathogenic",
                           "Uncertain_significance",
                           "Likely_benign","Benign")

patients <- vcf_pass %>% filter(GROUP == "Patient")
controls <- vcf_pass %>% filter(GROUP == "Control")

patients <- patients %>%
  mutate(GENE = sample(
    c(disease_genes, background_genes), nrow(patients),
    replace = TRUE,
    prob = c(rep(0.08, length(disease_genes)),
             rep(0.02, length(background_genes)))))

controls <- controls %>%
  mutate(GENE = sample(
    c(disease_genes, background_genes), nrow(controls),
    replace = TRUE,
    prob = c(rep(0.02, length(disease_genes)),
             rep(0.06, length(background_genes)))))

vcf_annotated <- bind_rows(patients, controls) %>%
  mutate(
    CONSEQUENCE = sample(consequence_types, nrow(.),
                         replace = TRUE,
                         prob = c(0.40,0.30,0.05,0.15,0.10)),
    CLIN_SIG    = sample(clinical_significance, nrow(.),
                         replace = TRUE,
                         prob = c(0.10,0.15,0.35,0.25,0.15)),
    POP_AF      = round(rbeta(nrow(.), 0.5, 5), 4)
  )

# SECTION 5: DISEASE ASSOCIATION
n_patients <- nrow(filter(vcf_annotated, GROUP == "Patient"))
n_controls <- nrow(filter(vcf_annotated, GROUP == "Control"))
results_list <- list()

for (gene in disease_genes) {
  pat_in  <- nrow(filter(vcf_annotated, GENE == gene, GROUP == "Patient"))
  con_in  <- nrow(filter(vcf_annotated, GENE == gene, GROUP == "Control"))
  pat_out <- n_patients - pat_in
  con_out <- n_controls - con_in
  tbl <- matrix(c(pat_in, con_in, pat_out, con_out), nrow = 2)
  ft  <- fisher.test(tbl)
  results_list[[gene]] <- tibble(
    GENE = gene, Patient = pat_in, Control = con_in,
    odds_ratio = round(ft$estimate, 3),
    p_value = round(ft$p.value, 4))
}

association_results <- bind_rows(results_list) %>%
  mutate(
    p_adjusted  = round(p.adjust(p_value, method = "bonferroni"), 4),
    significant = ifelse(p_adjusted < 0.05, "YES", "NO")
  ) %>%
  arrange(p_adjusted)

print(association_results)

# SECTION 6: SAVE ALL DATA
write.csv(vcf_data,
          paste0(project_folder, "/Project1_vcf_raw.csv"),
          row.names = FALSE)
write.csv(vcf_pass,
          paste0(project_folder, "/Project1_vcf_pass.csv"),
          row.names = FALSE)
write.csv(vcf_annotated,
          paste0(project_folder, "/Project1_vcf_annotated.csv"),
          row.names = FALSE)
write.csv(association_results,
          paste0(project_folder, "/Project1_association_results.csv"),
          row.names = FALSE)

cat("Project 1 Complete! All files saved.
")

