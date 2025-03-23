# sashimi command README

This document explains how to use the `ggsashimi.py` script to generate Sashimi plots for visualizing splicing events in RNA-seq data. The example command provided below demonstrates how to run the script with specific options and input files.

---

## Overview

**ggsashimi** is a Python-based tool that creates publication-quality Sashimi plots to visualize exon-exon junctions and read coverage over genomic intervals. This is particularly useful for analyzing alternative splicing patterns from RNA sequencing data.

---

## Prerequisites

Before running the script, ensure you have the following:

- **Python** installed on your system.
- Required Python libraries (e.g., `matplotlib`, `pandas`, `numpy`, etc.). Please refer to the documentation or `requirements.txt` file for the full list.
- Input files:
  - **BAM file**: Aligned and sorted RNA-seq reads (e.g., `SRR1259380.sorted.bam`).
  - **GTF file**: Genome annotation file (e.g., `dmel-all-r6.62.gtf` for *Drosophila melanogaster*).
  - **Palette file**: Custom color palette for plot styling (e.g., `palette.txt`).

---

## Command Breakdown

**Note** All the commands are to be run from the ./data directory.

Below is the example command and an explanation of its components:

```bash
python ../sashimi/ggsashimi.py \
  -b SRR1259380.sorted.bam \
  -c 2R:18757074-18761082 \
  -g ../database/dmel-all-r6.62.gtf \
  -o ../sashimi/sashimi \
  -M 10 \
  -C 3 \
  -O 0 \
  --shrink \
  --alpha 0.25 \
  --base-size=20 \
  --ann-height=5 \
  --height=5 \
  --width=18 \
  -P ../sashimi/palette.txt
```

### Parameter Explanations

- **`-b SRR1259380.sorted.bam`**  
  Specifies the input BAM file containing aligned RNA-seq reads.

- **`-c 2R:18757074-18761082`**  
  Sets the genomic coordinates to be visualized. Here, it represents a region on chromosome 2R from position 18757074 to 18761082.

- **`-g ../database/dmel-all-r6.62.gtf`**  
  Provides the genome annotation file (GTF format) necessary for identifying gene structures and exons.

- **`-o ../sashimi/sashimi`**  
  Specifies the output directory or file prefix for the generated plots.

- **`-M 10`**  
  Defines the maximum number of junctions or reads to be displayed (exact behavior might depend on the tool’s implementation).

- **`-C 3`**  
  Sets the minimum number of reads required to display a splice junction. Junctions with fewer reads will be omitted.

- **`-O 0`**  
  Specifies the offset value or other plotting parameter (e.g., control over how junction arcs are drawn).

- **`--shrink`**  
  Activates a mode that shrinks the plot to better fit the visual area (useful for long genomic regions).

- **`--alpha 0.25`**  
  Sets the transparency level of plotted elements to 0.25, which helps in overlaying multiple signals.

- **`--base-size=20`**  
  Determines the size of the base fonts in the plot for labels and annotations.

- **`--ann-height=5`**  
  Specifies the height allocated to the annotation track (gene models, exons, etc.).

- **`--height=5`**  
  Sets the overall height of the final plot.

- **`--width=18`**  
  Defines the overall width of the final plot.

- **`-P ../sashimi/palette.txt`**  
  Provides a custom palette file to define the colors used in the plot. This can be especially useful for maintaining a consistent visual style.

---

## Usage Example

To generate a Sashimi plot for a specific genomic region on chromosome 2R with the given settings, run the following command from your terminal:

```bash
python ../sashimi/ggsashimi.py -b SRR1259380.sorted.bam -c 2R:18757074-18761082 -g ../database/dmel-all-r6.62.gtf -o ../sashimi/sashimi -M 10 -C 3 -O 0 --shrink --alpha 0.25 --base-size=20 --ann-height=5 --height=5 --width=18 -P ../sashimi/palette.txt

# The above command output to PDF format and the below command output to SVG format
python ../sashimi/ggsashimi.py -b SRR1259380.sorted.bam -c 2R:18757074-18761082 -g ../database/dmel-all-r6.62.gtf -F svg -o ../sashimi/sashimi -M 10 -C 3 -O 0 --shrink --alpha 0.25 --base-size=20 --ann-height=5 --height=5 --width=18 -P ../sashimi/palette.txt
```

After execution, the output (plot image file(s)) will be saved in the designated output path (`../sashimi/sashimi`).

---

## Additional Notes

- **Customization:**  
  Adjust parameters such as `--alpha`, `--base-size`, `--ann-height`, `--height`, and `--width` to customize the appearance of the plot as needed.

  For different output formats, use `-F` followed by the desired format (e.g., `pdf`, `png`, `svg`, etc.).

- **Error Checking:**  
  Ensure that all input file paths are correct. If you encounter errors, check that the BAM file is properly sorted and indexed, and that the GTF and palette files are accessible.

- **Further Documentation:**  
  For more detailed information about each parameter or advanced usage, refer to the official ggsashimi documentation or use the help flag:
  ```bash
  python ../sashimi/ggsashimi.py --help
  ```