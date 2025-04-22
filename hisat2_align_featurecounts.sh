#!/bin/bash

# -----------------------------
# HISAT2 Alignment + featureCounts Pipeline
# -----------------------------

# HISAT2 index base name (without .ht2 extensions)
HISAT2_INDEX="Sus_genome_index"

# GTF annotation file
GTF_FILE="sus_scrofa.gtf"

# Number of threads to use
THREADS=8

# Output directory for alignments and counts
OUTPUT_DIR="./output"
mkdir -p "$OUTPUT_DIR"

# Loop through all *_1.fastq files (paired-end)
for FASTQ1 in *_1.fastq; do
    FASTQ2="${FASTQ1/_1.fastq/_2.fastq}"
    BASENAME=$(basename "$FASTQ1" _1.fastq)

    echo "Processing $BASENAME..."

    # Step 1: Align reads using HISAT2
    hisat2 -p $THREADS -x "$HISAT2_INDEX" -1 "$FASTQ1" -2 "$FASTQ2" -S "$OUTPUT_DIR/${BASENAME}.sam"

    # Step 2: Convert SAM to BAM
    samtools view -@ $THREADS -bS "$OUTPUT_DIR/${BASENAME}.sam" > "$OUTPUT_DIR/${BASENAME}.bam"

    # Step 3: Sort BAM file
    samtools sort -@ $THREADS -o "$OUTPUT_DIR/${BASENAME}_sorted.bam" "$OUTPUT_DIR/${BASENAME}.bam"

    # Step 4: Index BAM file
    samtools index "$OUTPUT_DIR/${BASENAME}_sorted.bam"

    # Step 5: Run featureCounts
    featureCounts -T $THREADS -a "$GTF_FILE" -o "$OUTPUT_DIR/${BASENAME}_counts.txt" -t exon -g gene_id "$OUTPUT_DIR/${BASENAME}_sorted.bam"

    # Clean up intermediate files
    rm "$OUTPUT_DIR/${BASENAME}.sam" "$OUTPUT_DIR/${BASENAME}.bam"

    echo "✔ Finished processing $BASENAME"
done

echo "✅ All samples processed."
