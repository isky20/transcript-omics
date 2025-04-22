This R script performs differential expression analysis using the edgeR package:
```
1- Reads a raw count matrix from a CSV file.
2- Extracts case and control sample columns based on user input.
3- Filters lowly expressed genes.
4- Normalizes the count data using TMM.
5- Builds a design matrix based on sample conditions.
6- Performs a likelihood ratio test for DE analysis.
7- Filters significant DEGs based on logFC and FDR thresholds.
8- Saves both normalized data and significant DEGs to CSV files.
```
```
Rscript run_edgeR_DEG.R \                          # Run the R script using Rscript
  --raw.reads.csv counts.csv \                     # Path to your input count matrix CSV file
  --colname.case Case1 \                           # Name of the first case sample column in the CSV
  --number.case.samples 3 \                        # Total number of case samples
  --named.case Case \                              # Label to assign to the case group
  --colname.control Control1 \                     # Name of the first control sample column in the CSV
  --number.control.samples 3 \                     # Total number of control samples
  --named.control Control \                        # Label to assign to the control group
  --number.logFC 1 \                               # Log2 fold-change threshold for filtering DEGs
  --number.FDR 0.05                                # FDR (adjusted p-value) threshold for filtering DEGs
```
