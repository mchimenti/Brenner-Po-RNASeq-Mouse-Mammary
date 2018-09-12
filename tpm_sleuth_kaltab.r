## Create TPM data from SALMON quantification for Brenner/Po project 
## Date: 09.12.2018
## Author: Michael Chimenti
## Organism: mm10
## Aligners: salmon
## Design: 
## Reps: 3

##########
## Imports
##########

#source("http://bioconductor.org/biocLite.R")
#biocLite("COMBINE-lab/wasabi")       #install wasabi the first time

#pseudoalignment analysis
library(wasabi)
library(sleuth)
library(biomaRt)

setwd("~/iihg/RNA_seq/brenner_lab/project_po_sept2018/")

#########################
## Salmon > HDF5 > Sleuth
#########################

base_dir <- "~/iihg/RNA_seq/brenner_lab/project_po_sept2018"

sal_dirs <- file.path(paste(base_dir, 
                            c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12"), "salmon", sep="/"))

prepare_fish_for_sleuth(sal_dirs)

## create a R matrix containing sample names and conditions from a text file
s2c <- read.table(file.path(base_dir, "samples.csv"), header = TRUE, stringsAsFactors=FALSE, sep = ',')

## add a column called "sal_dirs" containing paths to the data
s2c <- dplyr::mutate(s2c, path = sal_dirs)

## Get common gene names for transcripts

## this section queries Ensemble online database for gene names associated with transcript IDs
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")
t2g <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- rename(t2g, 'target_id' = 'ensembl_transcript_id', 'ens_gene' = 'ensembl_gene_id', 'ext_gene' = 'external_gene_name')

## Create the Sleuth Object

## create the sleuth object with the information matrix and the transcript to gene mapping
so <- sleuth_prep(s2c, ~ diet + treatment, target_mapping = t2g, aggregation_column = "ens_gene", extra_bootstrap_summary = FALSE)

kal_tab <- kallisto_table(so, include_covariates = FALSE)
#kal_tab <- kal_tab[, !names(kal_tab) %in% c("scaled_reads_per_base")]

## pivot long to wide
library(tidyr)
kal_tab_wide <- kal_tab %>% 
                as.tibble() %>% 
                select(c("target_id","sample","tpm")) %>% 
                spread(key = sample, value = tpm)

write.csv(file="kal_tab_wide_tpm.csv", kal_tab_wide)
