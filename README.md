This R script performs differential expression analysis 
using the edgeR package:
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
using the DEseq2 package:

```
1- Loads count matrix and sample metadata (study design).
2- Loops through each variable in the study design.
3- Runs DESeq2 for normalization and differential expression.
4- Filters significant genes (FDR ≤ 0.05).
5- Labels genes as up/downregulated.
6- Saves results for each comparison to output files.
```
Enrichment 
```
1- Load required libraries (dplyr, readxl, rbioapi).
2- Parse input arguments: CSV file name and STRING organism ID.
3- Read the input CSV containing gene symbols and regulation status.
4- Filter up-regulated genes and run STRING enrichment.
5- Save up-regulated enrichment results to a CSV file.
6- Filter down-regulated genes and run STRING enrichment.
7- Save down-regulated enrichment results to another CSV file.
8- If no results are found, save a message instead of a table.
```
