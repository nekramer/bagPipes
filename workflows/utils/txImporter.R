#!/usr/bin/R

## This script summarizes Salmon quantifications for RNA-Seq count matrix input to DESeq2

### INITIALIZE ===================================
## Load required libraries
library(tximeta)
library(readr)
library(tximport)

args <- commandArgs(trailingOnly=T)
samplesheet <- (args[1]) # samplesheet.csv
quantDir <- (args[2]) # output/quant/
runName <- (args[3]) # runName from snakemake
gtfFile <- args[4] # Optional gtf file to use rather than detected through tximeta

### READ IN ===================================
## Read in samplesheet file, name settings
coldata <- read_csv(samplesheet) |>
  dplyr::rename(names = sn) |> 
  # Paths to salmon files
  mutate(files = file.path(quantDir, names, "quant.sf"))

### RUN ===================================

# Transcript with inferred metadata or gtf 
if (!is.null(gtfFile)){

  ## Make TxDb if necessary
  txdb <- makeTxDbFromGFF(gtfFile, format = "gtf")
  k <- keys(txdb, keytype = "GENEID")
  df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")

  ## Assign proper TxDb columns to tx2gene object
  tx2gene <- df[, 2:1]
  se <- tximeta(coldata, tx2gene = tx2gene, skipMeta = TRUE, txOut = FALSE)
} else {
  # Import salmon transcript quants
  se <- tximeta(coldata)
}

# Convert to gene-level scaled transcripts
gse <- summarizeToGene(se)

# Write to file
save(gse, file = paste0(quantDir, "/", runName, "_gse.rda"))




