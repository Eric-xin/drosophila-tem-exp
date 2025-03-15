import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


def cis_regulation_F1(offspring_file, bias_file, imprint_indices, analysis_functions):
    """
    Reads F1 read counts, filters out biased and imprinting genes, and then for each temperature:
      - Subsets the data,
      - Constructs a design matrix based on allele (OregR vs. Samark),
      - Fits a negative binomial GLM per gene,
      - Extracts the allele effect as logFC and its p-value,
      - Adjusts the p-values to obtain FDR.
      
    Returns a DataFrame with cis-regulatory results per temperature.
    """
    # Read F1 read counts and bias info
    offspring = pd.read_csv(offspring_file, sep="\t", index_col="feature_id")
    gene_list = offspring.index.tolist()
    bias_genes = pd.read_csv(bias_file, sep="\t", index_col="feature_id")
    
    # Filter genes using the custom filtering function
    keep = offspring.apply(analysis_functions.filter_offspring_reads, axis=1)
    keep = keep & (~bias_genes["biased_sim"])
    for idx in imprint_indices:
        if idx in keep.index:
            keep.loc[idx] = False
    
    # Create grouping information for F1 samples.
    group = pd.DataFrame({
        "sample": [f"rep{i}" for i in list(range(1,13)) + list(range(13,25)) +
                             list(range(25,37)) + list(range(37,49))],
        "allele": np.tile(np.repeat(["OregR", "Samark"], 12), 2),
        "temperature": np.tile(np.repeat(['13','18','23','29'], 3), 4)
    })
    
    cis_results = {}
    # Define column subsets for each temperature (adjust these indices as needed)
    temp_cols = {
        '13': [1,2,3],
        '18': [4,5,6],
        '23': [7,8,9],
        '29': [9,10,11,19,20,21,33,34,35,45,46,47]
    }
    
    # For each temperature, fit GLMs gene-by-gene
    for temp, cols in temp_cols.items():
        subset = offspring.iloc[:, cols]
        group_sub = group.iloc[cols].reset_index(drop=True)
        # Build design matrix for allele effect (using one-hot encoding, with a constant)
        design = pd.get_dummies(group_sub["allele"], drop_first=True)
        design = sm.add_constant(design)
        
        gene_results = {}
        pvals = []
        gene_ids = []
        # Loop over each gene (row) in the subset
        for gene in subset.index:
            y = subset.loc[gene].values
            # Fit a negative binomial GLM to model counts ~ allele effect
            try:
                model = sm.GLM(y, design, family=sm.families.NegativeBinomial())
                result = model.fit()
                # Assume the allele effect is the second coefficient (after the constant)
                logFC = result.params[1]
                pval = result.pvalues[1]
            except Exception as e:
                logFC = np.nan
                pval = np.nan
            gene_results[gene] = {"logFC": logFC, "pval": pval}
            gene_ids.append(gene)
            pvals.append(pval)
        
        # Adjust p-values for multiple testing using FDR (Benjamini–Hochberg)
        _, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')
        # Create a DataFrame for this temperature with actual logFC and adjusted p-values (FDR)
        logFC_arr = [gene_results[gene]["logFC"] for gene in gene_ids]
        cis_table = pd.DataFrame({"logFC": logFC_arr, "FDR": pvals_adj}, index=gene_ids)
        cis_results[temp] = cis_table
    
    # Combine cis-regulatory results for all temperatures into one DataFrame
    combined = pd.DataFrame({"gene": gene_list})
    for temp in ['13', '18', '23', '29']:
        table = cis_results[temp]
        combined[f"logFC_t{temp}"] = np.nan
        combined[f"p_t{temp}"] = np.nan
        # Only assign results for genes that passed the filtering (keep == True)
        combined.loc[keep.index[keep], f"logFC_t{temp}"] = table["logFC"]
        combined.loc[keep.index[keep], f"p_t{temp}"] = table["FDR"]
    return combined

# --- Usage example ---
# Update these file paths to your actual file locations.
offspring_file = "./data/gene-readCounts-F1.txt"
bias_file = "./output/mapping-bias-genes-simulation.txt"
imprint_indices = [8021, 18699]

# Import your custom analysis functions module that at least provides filter_offspring_reads.
import analysis_functions

# Run cis-regulation analysis for F1
cis_F1 = cis_regulation_F1(offspring_file, bias_file, imprint_indices, analysis_functions)
# Save the results to a file
cis_F1.to_csv("./output/Cis-regulatory-F1-perTemperature.txt", sep="\t", index=False)