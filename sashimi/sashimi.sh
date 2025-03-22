#!/bin/bash
set -e  # Exit on error

# -----------------------------
# Configuration Variables
# -----------------------------
# Directories
DATADIR="data"           # For FASTQ, SAM, BAM, sorted BAM, output plots, etc.
DATABASEDIR="database"   # For reference genome, GTF, and HISAT2 index files

# Create directories if they do not exist
mkdir -p "$DATADIR" "$DATABASEDIR"

# Input FASTQ file (modify as needed)
FASTQ="${DATADIR}/sample.fastq"

# Reference genome (downloaded from FlyBase) and HISAT2 index variables
REFERENCE="${DATABASEDIR}/dmel-all-chromosome-r6.62.fasta"
REFERENCE_URL="https://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.62.fasta.gz"
REFERENCE_GZ="${REFERENCE}.gz"
HISAT2_INDEX="${DATABASEDIR}/dmel_r6.62_index"

# GTF annotation from FlyBase (for Drosophila melanogaster r6.62)
GTF_URL="https://ftp.flybase.org/genomes/dmel/current/gtf/dmel-all-r6.62.gtf.gz"
GTF_GZ="${DATABASEDIR}/dmel-all-r6.62.gtf.gz"
GTF="${DATABASEDIR}/dmel-all-r6.62.gtf"

# Region for sashimi plot (adjust chromosome and coordinates as needed)
REGION="chr12:56541490-56572047"

# Intermediate and output file names (located in DATADIR)
SAM_FILE="${DATADIR}/sample.sam"
BAM_FILE="${DATADIR}/sample.bam"
SORTED_BAM="${DATADIR}/sample.sorted.bam"
OUT_PLOT="${DATADIR}/sample_sashimi_plot.png"

# -----------------------------
# Step 1: Download and Unzip the Reference Genome
# -----------------------------
if [ ! -f "${REFERENCE}" ]; then
    echo "Reference genome not found. Downloading and unzipping..."
    wget "$REFERENCE_URL" -O "$REFERENCE_GZ"
    gunzip -k "$REFERENCE_GZ"  # Keep the gzipped file for reference if desired
fi

# -----------------------------
# Step 2: Build the HISAT2 Index (if not already built)
# -----------------------------
if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
    echo "Building HISAT2 index from ${REFERENCE}..."
    hisat2-build "${REFERENCE}" "${HISAT2_INDEX}"
fi

# -----------------------------
# Step 3: Quality Control on FASTQ
# -----------------------------
echo "Running FastQC on ${FASTQ}..."
fastqc "${FASTQ}"

# -----------------------------
# Step 4: Align Reads Using HISAT2
# -----------------------------
echo "Aligning reads with HISAT2..."
hisat2 -p 8 -x "${HISAT2_INDEX}" -U "${FASTQ}" -S "${SAM_FILE}"

# -----------------------------
# Step 5: Convert SAM to BAM
# -----------------------------
echo "Converting SAM to BAM..."
samtools view -bS "${SAM_FILE}" > "${BAM_FILE}"

# -----------------------------
# Step 6: Sort BAM File
# -----------------------------
echo "Sorting BAM file..."
samtools sort "${BAM_FILE}" -o "${SORTED_BAM}"

# -----------------------------
# Step 7: Index Sorted BAM File
# -----------------------------
echo "Indexing sorted BAM file..."
samtools index "${SORTED_BAM}"

# -----------------------------
# Step 8: Download and Prepare GTF Annotation
# -----------------------------
if [ ! -f "${GTF}" ]; then
    if [ ! -f "${GTF_GZ}" ]; then
        echo "Downloading GTF annotation from FlyBase..."
        wget "$GTF_URL" -O "${GTF_GZ}"
    fi
    echo "Decompressing GTF annotation..."
    gunzip -k "${GTF_GZ}"  # Keep the original gzipped file
fi

# -----------------------------
# Step 9: Generate Sashimi Plot with ggsashimi
# -----------------------------
echo "Generating sashimi plot..."
ggsashimi -b "${SORTED_BAM}" -c "${REGION}" -g "${GTF}" -o "${OUT_PLOT}"

echo "Sashimi plot generated: ${OUT_PLOT}"