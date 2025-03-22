# HISAT2 Alignment Workflow for Drosophila Sashimi Plot

---

## Step 1: Download and Unzip the Reference Genome

Fetch the Drosophila melanogaster reference genome from [FlyBase](https://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.62.fasta.gz) by running:

```sh
wget https://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.62.fasta.gz
gunzip dmel-all-chromosome-r6.62.fasta.gz
```

---

## Step 2: Build the HISAT2 Index

Create a HISAT2 index from the downloaded reference genome using:

```sh
hisat2-build dmel-all-chromosome-r6.62.fasta dmel_r6.62_index
```

---

## Step 3: Align Your FASTQ File Using HISAT2

Align your FASTQ file against the reference genome with the generated HISAT2 index:

```sh
hisat2 -p 8 -x dmel_r6.62_index -U sample.fastq -S sample.sam
```

For example, if your files are located in specific directories, you might run:

```sh
hisat2 -p 8 -x ../database/index/dmel_r6.62_index -U SRR1259380_1.fastq -S SRR1234567.sam
```

---

### Example Output

```
31186966 reads; of these:
  31186966 (100.00%) were unpaired; of these:
    927776 (2.97%) aligned 0 times
    29657485 (95.10%) aligned exactly 1 time
    601705 (1.93%) aligned >1 times
97.03% overall alignment rate
```