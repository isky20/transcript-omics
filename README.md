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
