Below is a revised version of the workflow instructions and accompanying notes:

---

## Sashimi Workflow

1. **Quality Control:**  
   Run FastQC on your FASTQ file:  
   ```bash
   fastqc sample.fastq
   ```

2. **Build HISAT2 Index:**  
   Create a HISAT2 index from your reference genome:  
   ```bash
   hisat2-build reference.fa genome
   ```

3. **Align Reads:**  
   Align your reads using HISAT2:  
   ```bash
   hisat2 -p 8 -x genome -U sample.fastq -S sample.sam
   ```

4. **Convert SAM to BAM:**  
   Convert the SAM file to BAM format:  
   ```bash
   samtools view -bS sample.sam > sample.bam
   ```

5. **Sort BAM File:**  
   Sort the BAM file:  
   ```bash
   samtools sort sample.bam -o sample.sorted.bam
   ```

6. **Index BAM File:**  
   Index the sorted BAM file:  
   ```bash
   samtools index sample.sorted.bam
   ```

7. **Download GTF Annotation:**  
   Download the annotation file from FlyBase:  
   [FlyBase dmel-all-r6.62.gtf.gz](https://ftp.flybase.org/genomes/dmel/current/gtf/dmel-all-r6.62.gtf.gz)

8. **Generate Sashimi Plot:**  
   Use ggsashimi to create the plot for a specific region:  
   ```bash
   ggsashimi -b sample.subset.bam -c chr12:56541490-56572047 -g annotation.gtf -o sample_sashimi_plot.png
   ```

---

## Additional Notes
- **Running ggsashimi on Different Regions:**  
  Example commands for different genomic regions:  
  ```bash
  python ../sashimi/ggsashimi.py -b SRR1259380.sorted.bam -c 3R:20669571-20672821 -g ../database/dmel-all-r6.62.gtf -o ../sashimi/sashimi
  ```
  ```bash
  python ../sashimi/ggsashimi.py -b SRR1259380.sorted.bam -c 2R:14600000-14650000 -g ../database/dmel-all-r6.62.gtf -o ../sashimi/sashimi
  ```

- **Final Figure Generation:**  
  Create a final figure using a dedicated script:  
  ```bash
  python ../sashimi/figure.py -b SRR1259380.sorted.bam -c 2R:14600000-14650000 -o ../sashimi/sashimi_plot.png
  ```