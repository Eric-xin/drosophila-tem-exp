# Download Sequencing Data from the NCBI Sequence Read Archive (SRA) Using the SRA Toolkit

1. **Install the SRA Toolkit**

2. **Configure the Toolkit (Optional)**

    **Initial Setup:**
    After installation, you might need to run the configuration to set the default download directory. For example:

    ```sh
    vdb-config --interactive
    ```

    This will launch an interactive configuration where you can set your preferred settings.

---

3. **Download Data Using an Accession Number**

    - **Using prefetch:**
      The `prefetch` command downloads the SRA file associated with a given accession number. For example, to download the SRA file for accession `SRR1259380`:

      ```sh
      prefetch SRR1259380
      ```

      This command downloads the file to the default SRA download directory.

---

4. **Convert the SRA File to FASTQ Format**

    - **Using fasterq-dump:**
      Once the SRA file is downloaded, you can convert it to FASTQ format using `fasterq-dump`. For example:

      ```sh
      fasterq-dump SRR1259380
      ```

      This command converts the SRA file into one or more FASTQ files in your current working directory.

    - **Paired-End Data:**
      If your data is paired-end, you can use the `--split-files` option to split the reads into separate files:

      ```sh
      fasterq-dump SRR1259380 --split-files
      ```

---

5. **Verify and Use Your FASTQ Files**

    After the conversion, you should see one or more FASTQ files in your directory. These files can now be used for downstream analysis such as quality control, alignment, or sashimi plotting.

---

## Summary

1. Install the SRA Toolkit.
2. Configure the toolkit using `vdb-config` if needed.
3. Download your SRA file with `prefetch <accession>`.
4. Convert the SRA file to FASTQ with `fasterq-dump <accession>` (using `--split-files` for paired-end reads).