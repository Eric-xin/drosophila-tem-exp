# Downloading and Converting NCBI|SRR Data to Sorted BAM

This guide will walk you through downloading SRR1234567 sequencing data using the SRA Toolkit, converting it to FASTQ, aligning the reads with HISAT2, and converting the alignment to a sorted BAM file using Samtools.

## Prerequisites

Ensure you have the following installed and accessible in your PATH:

- [SRA Toolkit](https://github.com/ncbi/sra-tools) (includes `prefetch` and `fasterq-dump`)
- [HISAT2](https://daehwankimlab.github.io/hisat2/)
- [Samtools](http://www.htslib.org/)

## Step 1: Download the Data with SRA Toolkit

First, use `prefetch` to download the SRR1234567 dataset from the SRA database:

```bash
prefetch SRR1234567
```

The downloaded file will be placed in the SRA Toolkit’s default directory (often `~/ncbi/public/sra/`).

## Note for Mac Users
use this command to avoid malware issue
```sh
sudo xattr -r -d com.apple.quarantine /Users/ericxin/code/sratoolkit/bin/*
```

## Step 2: Convert SRA to FASTQ

Next, convert the downloaded SRA file to FASTQ format using `fasterq-dump` (which is generally faster than `fastq-dump`):

```bash
fasterq-dump SRR1234567
```

This command will create one or more FASTQ files (for single-end or paired-end data) in your current directory. For paired-end data, you might see files like:
- `SRR1234567_1.fastq`
- `SRR1234567_2.fastq`

For single-end data, you’ll have a single `SRR1234567.fastq` file.

## Step 3: Align Reads with HISAT2

Use HISAT2 to align your FASTQ file(s) to your reference genome. For example, if you are working with single-end data and have built your HISAT2 index (see [HISAT2 Documentation](https://daehwankimlab.github.io/hisat2/manual/)), run:

```bash
hisat2 -p 8 -x /path/to/hisat2/index -U SRR1234567.fastq -S SRR1234567.sam
```

For paired-end data, use:

```bash
hisat2 -p 8 -x /path/to/hisat2/index -1 SRR1234567_1.fastq -2 SRR1234567_2.fastq -S SRR1234567.sam
```

This command creates an alignment file in SAM format.

## Step 4: Convert SAM to Sorted BAM

Convert the SAM file to BAM format, sort the BAM file, and create an index using Samtools.

1. **Convert SAM to BAM:**

   ```bash
   samtools view -bS SRR1234567.sam > SRR1234567.bam
   ```

2. **Sort the BAM File:**

   ```bash
   samtools sort SRR1234567.bam -o SRR1234567.sorted.bam
   ```

3. **Index the Sorted BAM File (Optional, but useful for visualization):**

   ```bash
   samtools index SRR1234567.sorted.bam
   ```

## Final Notes

- **Quality Control:**  
  You may wish to run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on the FASTQ files before alignment to assess quality.

- **File Organization:**  
  For convenience, you can organize your files into directories (e.g., `data/` for FASTQ, SAM, and BAM files and `database/` for reference and index files).

- **Parameters:**  
  Adjust the number of threads (`-p 8` in HISAT2) and file paths as needed based on your system and data.

By following these steps, you'll have a sorted BAM file (`SRR1234567.sorted.bam`) that you can use for downstream analysis such as generating sashimi plots or other genomic visualizations.