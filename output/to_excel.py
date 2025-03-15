import pandas as pd
import numpy as np

"""
This script reads a TSV file containing processed allele counts, sorts the data by chromosome arm and start position, and rounds the numeric values to three significant figures.
The rounded data is then written to an Excel file for further analysis.
"""

def round_to_sig(x, sig=3):
    """
    Rounds a number x to a given number of significant figures.
    If x is not a number, returns it unchanged.
    """
    try:
        # Convert to float then format to three significant figures and convert back to float.
        return float(f"{float(x):.{sig}g}")
    except (ValueError, TypeError):
        return x

df = pd.read_csv("./output/processed_allele_counts.tsv", delimiter="\t")  # or delimiter="," if it's comma separated

df_sorted = df.sort_values(by=["arm", "start"])

# Apply the rounding function to numeric columns
numeric_cols = df_sorted.select_dtypes(include=[np.number]).columns
for col in numeric_cols:
    df_sorted[col] = df_sorted[col].apply(round_to_sig)

# Write the modified DataFrame to an Excel file.
df_sorted.to_excel("./output/processed_allele_counts.tsv.xlsx", index=False)