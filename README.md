# 🧬 RNA-Seq Read Alignment, Quantification, and Enrichment Pipeline

This pipeline is a simplified RNA-Seq analysis workflow that includes read alignment, gene quantification, differential expression analysis, and functional enrichment analysis.

---

## 📂 Input Requirements
- Paired-end FASTQ files (`*_R1.fastq*`, `*_R2.fastq*`)
- Reference genome STAR index
- GTF gene annotation file
- Study design file (CSV, describing sample groups)
- STRING taxonomy ID (e.g., `9606` for human)

---

## 🔹 Pipeline Steps

### 📁 1. Prepare Output Directory
- Create a directory to save all output files.

### 🔎 2. Find FASTQ Files
- Search and identify paired-end FASTQ files based on standard naming (`*_R1.fastq*` and `*_R2.fastq*`).

### 🧬 3. Align Reads Using STAR
- Align reads to the reference genome using **STAR**.
- **Output**: Sorted BAM files for each sample.

### 📊 4. Generate Gene Counts with featureCounts
- Use **featureCounts** to quantify reads mapped to genes based on the GTF annotation file.
- **Output**: A combined gene counts matrix (`all_samples_gene_counts.txt`).

---

# 🧬 Differential Expression Analysis with DESeq2

This stage identifies differentially expressed genes (DEGs) between conditions based on normalized gene counts.

---

## 📂 Input Requirements
- Gene counts matrix (`count_p10_mrna.csv`)
- Sample metadata/study design file (`sd.csv`)

---

## 🔹 Pipeline Steps

### 📁 1. Load Gene Counts and Study Design
- Load raw gene counts and corresponding sample group information.

### 🔎 2. Loop Over Each Condition
- For each condition/factor in the study design, create a model for differential expression analysis.

### 📈 3. Normalize and Perform DESeq2 Analysis
- Normalize library sizes and perform differential expression testing.

### 🧹 4. Filter Significant DEGs
- Extract genes with adjusted p-value ≤ 0.05.

### 💾 5. Save DESeq2 Results
- Save significant DEGs into `significant_results_*.csv`.
- If no DEGs are detected, create a `no_significant_results_*.csv`.

---

# 🧬 Functional Enrichment Analysis with STRING

This stage identifies enriched biological processes, pathways, and functions among upregulated and downregulated genes.

---

## 📂 Input Requirements
- DESeq2 results file (`significant_results_*.csv`)
- STRING organism Taxonomy ID (e.g., `9606` for human, `9823` for pig)

---

## 🔹 Pipeline Steps

### 📁 1. Load DEG Results
- Read in the DEGs from DESeq2 output.

### ➗ 2. Separate Upregulated and Downregulated Genes
- Identify "upregulated" and "downregulated" genes based on the `regulation` column.

### 🔥 3. Perform STRING Enrichment Analysis
- Perform functional enrichment analysis separately for:
  - Upregulated genes
  - Downregulated genes
- Using the **STRING database** through **rbioapi**.

### 💾 4. Save Enrichment Results
- Save enrichment results into:
  - `*_up.csv` (upregulated enrichment)
  - `*_down.csv` (downregulated enrichment)

---

# 🚀 Example Commands

### 🔹 Run the Python STAR + featureCounts Pipeline

```
python3 run_rna_pipeline.py \
  --fastq-dir test_data \
  --index-dir star_index \
  --gtf-file annotation.gtf \
  --output-dir output \
  --threads 2
```

```
Rscript run_edgeR_DEG.R \
  --raw.reads.csv counts.csv \
  --colname.case Case1 \
  --number.case.samples 3 \
  --named.case Case \
  --colname.control Control1 \
  --number.control.samples 3 \
  --named.control Control \
  --number.logFC 1 \
  --number.FDR 0.05
```
