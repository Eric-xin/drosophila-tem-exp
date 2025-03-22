#!/usr/bin/env python3
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import argparse

def get_coverage(bamfile, chrom, start, end):
    """Return a numpy array with coverage from start to end on chrom."""
    # Initialize coverage vector with zeros
    coverage = np.zeros(end - start + 1, dtype=int)
    for pileupcolumn in bamfile.pileup(chrom, start, end, truncate=True):
        pos = pileupcolumn.reference_pos
        if start <= pos <= end:
            coverage[pos - start] = pileupcolumn.nsegments
    return coverage

def get_junctions(bamfile, chrom, start, end):
    """
    Parse through alignments in the region and extract splice junctions.
    Returns a dictionary with keys as (junction_start, junction_end) and values as counts.
    """
    junctions = {}
    # Iterate over reads overlapping the region
    for read in bamfile.fetch(chrom, start, end):
        # Only consider primary alignments
        if read.is_secondary or read.is_supplementary:
            continue
        # Check for spliced reads: an 'N' operation in the CIGAR string indicates an intron
        if read.cigartuples is None:
            continue
        ref_pos = read.reference_start
        for (op, length) in read.cigartuples:
            # CIGAR op 3 corresponds to 'N' (skipped region, i.e. intron)
            if op == 3:
                junction_start = ref_pos
                junction_end = ref_pos + length
                key = (junction_start, junction_end)
                junctions[key] = junctions.get(key, 0) + 1
            # For operations that consume reference, update position
            if op in (0, 2, 3, 7, 8):
                ref_pos += length
    return junctions

def plot_sashimi(chrom, start, end, coverage, junctions, out_file):
    """Create a sashimi-like plot from coverage and junctions data."""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot read coverage
    positions = np.arange(start, end + 1)
    ax.fill_between(positions, coverage, color='lightgray', step='mid', label='Coverage')
    ax.set_xlabel(f"{chrom} position")
    ax.set_ylabel("Read Coverage")
    ax.set_title(f"Sashimi Plot: {chrom}:{start}-{end}")

    # Plot junction arcs
    # Set parameters for arcs
    arc_color = 'dodgerblue'
    for (jstart, jend), count in junctions.items():
        # Only plot junctions that are fully inside the region
        if jstart < start or jend > end:
            continue
        mid = (jstart + jend) / 2
        width = jend - jstart
        # Set arc height proportional to the log of count, or simply proportional to count.
        arc_height = 0.1 * count  # adjust factor as needed

        # Create an arc: center at (mid, current coverage at mid), width equals jend - jstart, height as defined.
        # We use a FancyArrowPatch to simulate an arc.
        # Alternatively, use matplotlib.patches.Arc if you want a simple arc shape.
        arc = patches.Arc((mid, max(coverage)), width=width, height=arc_height * 10,
                          angle=0, theta1=0, theta2=180, color=arc_color, lw=2)
        ax.add_patch(arc)
        # Add a text label at the top center of the arc with count
        ax.text(mid, max(coverage) + arc_height * 10 + 1, str(count),
                ha='center', va='bottom', color=arc_color, fontsize=10)

    ax.legend()
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    print(f"Plot saved as {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Custom Sashimi Plot from BAM and GTF.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-c", "--coordinates", required=True,
                        help="Coordinates in the format chrom:start-end (e.g., 2R:14600000-14650000)")
    parser.add_argument("-o", "--out", required=True, help="Output image file (e.g., plot.png)")
    args = parser.parse_args()

    # Parse coordinates
    try:
        chrom, pos = args.coordinates.split(":")
        start, end = map(int, pos.split("-"))
    except Exception as e:
        parser.error("Coordinates must be in the format chrom:start-end")
    
    # Open BAM file
    bamfile = pysam.AlignmentFile(args.bam, "rb")

    # Get coverage vector for the region
    coverage = get_coverage(bamfile, chrom, start, end)
    # Get splice junctions from reads in the region
    junctions = get_junctions(bamfile, chrom, start, end)

    bamfile.close()

    # Plot and save the sashimi-like plot
    plot_sashimi(chrom, start, end, coverage, junctions, args.out)

if __name__ == "__main__":
    main()