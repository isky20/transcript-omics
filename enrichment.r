#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)      # For data manipulation
library(readxl)     # For reading Excel files (not used here but loaded in case needed)
library(rbioapi)    # For STRING enrichment analysis using STRING API

# Define the enrichment analysis function
perform_enrichment_analysis <- function(file_name, organism_id) {
  
  # Check if the input file exists
  if (!file.exists(file_name)) {
    stop("File does not exist: ", file_name)
  }
  
  # Read CSV file containing gene symbols and regulation info
  df <- read.csv(file_name)
  
  # ---- Enrichment analysis for up-regulated genes ----
  
  # Extract up-regulated gene symbols
  up_list <- df$Symbol[grepl("up", df$regulation, ignore.case = TRUE)]
  
  # Perform STRING enrichment analysis using rbioapi
  enrichmentup <- rbioapi::rba_string_enrichment(up_list, organism_id, split_df = FALSE)
  
  # Check if results are found; if so, save to CSV
  if (nrow(enrichmentup) == 0) {
    cat("No enrichment found for up-regulated genes in", file_name, "\n")
  } else {
    subset_dfUP <- enrichmentup[, c("term", "description", "category", "number_of_genes")]
    write.csv(subset_dfUP, paste0(gsub(".csv", "", file_name), "_up.csv"), row.names = FALSE)
  }
  
  # ---- Enrichment analysis for down-regulated genes ----
  
  # Extract down-regulated gene symbols
  down_list <- df$Symbol[grepl("down", df$regulation, ignore.case = TRUE)]
  
  # Perform STRING enrichment analysis using rbioapi
  enrichment <- rbioapi::rba_string_enrichment(down_list, organism_id, split_df = FALSE)
  
  # Check if results are found; if so, save to CSV
  if (nrow(enrichment) == 0) {
    cat("No enrichment found for down-regulated genes in", file_name, "\n")
  } else {
    subset_df <- enrichment[, c("term", "description", "category", "number_of_genes")]
    write.csv(subset_df, paste0(gsub(".csv", "", file_name), "_down.csv"), row.names = FALSE)
  }
}

# ---- Command-line execution block ----

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Require two arguments: input file name and organism ID
if (length(args) < 2) {
  stop("Usage: Rscript perform_enrichment_analysis.R <file_name> <organism_id>")
}

# Assign arguments to variables
file_name <- args[1]               # Path to input CSV file
organism_id <- as.numeric(args[2]) # STRING taxonomy ID for the organism

# Run the function with provided arguments
perform_enrichment_analysis(file_name, organism_id)
