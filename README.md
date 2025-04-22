This R script performs differential expression analysis 

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
RNA pipline with STAR
```
Parse arguments – Get paths and settings from the user.
Detect FASTQ files – Identify _R1 and _R2 pairs (or single-end).
Run STAR – Align each sample to the reference genome and generate sorted BAM files.
Collect BAMs – Store output BAM files for downstream use.
Run featureCounts – Quantify gene expression from BAMs using the provided GTF.
Save Output – Store gene count matrix in the specified output directory.
```
```
python3 run_rna_pipeline.py \
  --fastq-dir test_data \
  --index-dir star_index \
  --gtf-file annotation.gtf \
  --output-dir output \
  --threads 2
```
