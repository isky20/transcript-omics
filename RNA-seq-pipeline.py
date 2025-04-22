#!/usr/bin/env python3

import os
import subprocess
import argparse
import glob

def run_star(fastq1, fastq2, index_dir, output_dir, threads, is_paired, sample_name):
    """Run STAR alignment for one sample."""

    # Create output directory for the sample
    sample_output_dir = os.path.join(output_dir, sample_name)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Build STAR command
    star_cmd = [
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", index_dir,
        "--outFileNamePrefix", os.path.join(sample_output_dir, "STAR_"),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--quantMode", "GeneCounts"
    ]

    # Add FASTQ input (single or paired)
    if is_paired:
        star_cmd += ["--readFilesIn", fastq1, fastq2]
    else:
        star_cmd += ["--readFilesIn", fastq1]

    # Run STAR
    print(f"Running STAR for {sample_name}...")
    subprocess.run(star_cmd, check=True)

    # Return path to sorted BAM
    return os.path.join(sample_output_dir, "STAR_Aligned.sortedByCoord.out.bam")


def run_featurecounts(bam_files, gtf_file, output_file, threads, is_paired):
    """Run featureCounts on a list of BAM files."""

    # Build featureCounts command
    featurecounts_cmd = [
        "featureCounts",
        "-T", str(threads),
        "-a", gtf_file,
        "-o", output_file,
        "-t", "exon",
        "-g", "gene_id"
    ]

    if is_paired:
        featurecounts_cmd.append("-p")

    featurecounts_cmd += bam_files

    # Run featureCounts
    print("Running featureCounts on all samples...")
    subprocess.run(featurecounts_cmd, check=True)


def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="RNA-Seq pipeline using STAR and featureCounts.")
    parser.add_argument("--fastq-dir", required=True, help="Directory containing FASTQ files.")
    parser.add_argument("--index-dir", required=True, help="STAR genome index directory.")
    parser.add_argument("--gtf-file", required=True, help="GTF annotation file.")
    parser.add_argument("--output-dir", required=True, help="Directory for output files.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads.")
    args = parser.parse_args()

    # Detect FASTQ files
    fastq1_files = sorted(glob.glob(os.path.join(args.fastq_dir, "*_R1*.fastq*")))
    fastq2_files = sorted(glob.glob(os.path.join(args.fastq_dir, "*_R2*.fastq*")))
    bam_files = []

    # Process each FASTQ sample
    for fastq1 in fastq1_files:
        sample_name = os.path.basename(fastq1).split("_R1")[0]
        fastq2 = fastq1.replace("_R1", "_R2")
        is_paired = os.path.exists(fastq2) if fastq2_files else False

        # Run STAR and collect BAM file
        bam_file = run_star(fastq1, fastq2, args.index_dir, args.output_dir, args.threads, is_paired, sample_name)
        bam_files.append(bam_file)

    # Run featureCounts on all aligned BAMs
    featurecounts_output = os.path.join(args.output_dir, "all_samples_gene_counts.txt")
    run_featurecounts(bam_files, args.gtf_file, featurecounts_output, args.threads, is_paired)


if __name__ == "__main__":
    main()
