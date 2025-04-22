# Load required libraries
library(edgeR)      # For differential expression analysis
library(dplyr)      # For data manipulation
library(argparse)   # For parsing command-line arguments

# Function to get Differentially Expressed Genes (DEGs)
get_DEG <- function(raw.reads.csv, colname.case, number.case.samples, named.case,
                    colname.control, number.control.samples, named.control, number.logFC, number.FDR) {
  
  # Load count data from CSV file
  rawcounts <- read.csv(raw.reads.csv, header = TRUE, row.names = 1)
  
  # Find index of case and control samples based on column names
  case_start <- which(colnames(rawcounts) == colname.case)
  case_end <- case_start + number.case.samples - 1
  control_start <- which(colnames(rawcounts) == colname.control)
  control_end <- control_start + number.control.samples - 1
  
  # Subset counts to include only case and control columns
  counts <- rawcounts[, c(case_start:case_end, control_start:control_end)]
  
  # Create a DGEList object for edgeR analysis
  dge <- DGEList(counts = counts)
  
  # Filter out lowly expressed genes (genes expressed in at least half the samples)
  keep <- rowSums(cpm(dge) > 0) >= ceiling(ncol(counts) / 2)
  dge <- dge[keep, ]
  
  # Normalize counts using TMM method
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Export normalized CPM data to CSV
  norm_data <- cpm(dge, normalized.lib.sizes = TRUE)
  write.csv(norm_data, file = paste("normal_data_", named.case, "_", named.control, ".csv"), quote = TRUE, row.names = TRUE)
  
  # Create metadata (sample labels for case and control groups)
  metadata <- data.frame(
    sample_id = colnames(counts),
    condition = c(rep(named.case, times = number.case.samples), rep(named.control, times = number.control.samples))
  )
  metadata$condition <- relevel(factor(metadata$condition), ref = named.control)  # Set reference level
  
  # Design matrix for modeling
  design <- model.matrix(~ condition, data = metadata)
  
  # Estimate dispersion for the count data
  dge <- estimateDisp(dge, design)
  
  # Fit GLM model and perform likelihood ratio test
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit)
  
  # Extract all DE results
  de_genes <- topTags(lrt, n = 99999999)$table
  
  # Keep only genes with p-value < 0.05, sort by logFC
  de_genes_sig <- arrange(de_genes[de_genes$PValue < 0.05, ], logFC)
  
  # Annotate results with regulation direction (up/down)
  tibble_data <- de_genes_sig %>%
    mutate(regulation = ifelse(logFC > 0, paste("up in", named.case), paste("down in", named.case)))
  
  # Convert to data frame and filter based on thresholds
  result <- as.data.frame(tibble_data)
  filtered <- subset(result, abs(logFC) > number.logFC & FDR < number.FDR)
  
  # Save significant DEGs to file
  write.csv(filtered, file = paste(named.case, named.control, "logFC", number.logFC, "FDR", number.FDR, "DEGsNEW.csv", sep = "_"), quote = TRUE, row.names = TRUE)
}

# Wrapper function to run from command line
run_with_args <- function() {
  parser <- ArgumentParser(description = "Differential Expression Analysis using edgeR")
  
  # Define expected arguments
  parser$add_argument("--raw.reads.csv", required = TRUE, help = "Path to the raw counts CSV file")
  parser$add_argument("--colname.case", required = TRUE, help = "Column name for the first case sample")
  parser$add_argument("--number.case.samples", type = "integer", required = TRUE, help = "Number of case samples")
  parser$add_argument("--named.case", required = TRUE, help = "Name for the case group")
  parser$add_argument("--colname.control", required = TRUE, help = "Column name for the first control sample")
  parser$add_argument("--number.control.samples", type = "integer", required = TRUE, help = "Number of control samples")
  parser$add_argument("--named.control", required = TRUE, help = "Name for the control group")
  parser$add_argument("--number.logFC", type = "numeric", required = TRUE, help = "logFC threshold")
  parser$add_argument("--number.FDR", type = "numeric", required = TRUE, help = "FDR threshold")
  
  # Parse arguments from the command line
  args <- parser$parse_args()
  
  # Call main DEG function with parsed inputs
  get_DEG(args$raw.reads.csv, args$colname.case, args$number.case.samples, args$named.case,
          args$colname.control, args$number.control.samples, args$named.control,
          args$number.logFC, args$number.FDR)
}

# Run the script when executed
run_with_args()
