# Load necessary libraries
library(DESeq2)   # For performing differential expression analysis
library(tidyr)    # For data manipulation

# Get command-line arguments (not used in current script, kept for future flexibility)
args <- commandArgs(trailingOnly = TRUE)

# Define file paths manually
input.count.matrix <- "count_p10_mrna.csv"   # CSV file with raw count data (genes x samples)
input.study.design <- "sd.csv"               # CSV file with sample metadata (study design)
output.directory <- "out"                    # Directory to save the output results

# Load count matrix; row names = gene names
data.matrix <- read.csv(file = input.count.matrix, header = TRUE, sep = ",", row.names = "Name", stringsAsFactors = FALSE)

# Load study design matrix; row names = sample names
study.design.matrix <- read.csv(input.study.design, header = TRUE, row.names = "samplename")

# Loop over each factor (column) in the study design matrix
for (i in colnames(study.design.matrix)) {
  
  formula <- as.formula(paste("~", i))  # Create a model formula dynamically using the current column
  
  # Create a DESeq2 dataset object using counts and design
  DE.object <- DESeqDataSetFromMatrix(
    countData = data.matrix,
    colData = study.design.matrix,
    design = formula
  )
  
  # Estimate size factors for library normalization
  DE.object <- estimateSizeFactors(DE.object)
  
  # Perform DESeq2 analysis
  DE.object <- DESeq(DE.object)
  
  # Extract results for the current factor comparison
  results <- results(DE.object)
  
  # Identify the reference and comparison groups (factor levels)
  opposite_level <- levels(colData(DE.object)[[i]])[2]
  ref_level <- levels(colData(DE.object)[[i]])[1]
  
  # Filter significant results (adjusted p-value ≤ 0.05)
  res_filtered <- results[!is.na(results$padj) & results$padj <= 0.05, ]
  
  # Check if any significant genes were found
  if (nrow(res_filtered) == 0) {
    # If no DEGs found, write a message to file
    message <- paste("No significant results between", opposite_level, "and", ref_level)
    write.csv(data.frame(Message = message), 
              file.path(output.directory, paste0("no_significant_results_", i, ".csv")))
  } else {
    # Add a column indicating regulation direction (up/down)
    res_filtered$regulation <- ifelse(
      res_filtered$log2FoldChange > 0,
      paste("upregulated", opposite_level, "Vs.", ref_level),
      paste("downregulated", opposite_level, "Vs.", ref_level)
    )
    
    # Save the filtered results to CSV
    write.csv(res_filtered, 
              file.path(output.directory, paste0("significant_results_", i, "_", opposite_level, "Vs", ref_level, ".csv")))
  }
}
