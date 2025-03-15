#!/usr/bin/env python
import pandas as pd
import numpy as np

# Import your custom analysis functions module.
# This module must implement the function calc_sim_factor(row)
import analysis_functions

def simulation_correction(simulation_file, bias_output_file, analysis_functions):
    """
    Reads simulation data from a file containing columns:
    feature_id, Parent1, and Parent2. Calculates a simulation factor per gene using
    calc_sim_factor, flags genes as biased if the factor is not equal to 1 or is missing,
    and writes the bias information to an output file.
    
    Parameters:
      simulation_file (str): Path to the simulation file.
      bias_output_file (str): Path to write the bias information.
      analysis_functions: Module containing the custom calc_sim_factor function.
      
    Returns:
      simFactors (pd.Series): Series of simulation factors per gene.
      bias_df (pd.DataFrame): DataFrame with columns 'feature_id' and 'biased_sim'.
    """
    # Read the simulation file.
    sim = pd.read_csv(simulation_file, sep="\t")
    sim.set_index("feature_id", inplace=True)
    
    # Calculate simulation factors row-by-row using the custom function.
    simFactors = sim.apply(lambda row: analysis_functions.calc_sim_factor(row), axis=1)
    
    # Flag genes with a simulation factor not equal to 1 or that are missing.
    biased_sim = (simFactors != 1) | (simFactors.isna())
    bias_df = pd.DataFrame({
        "feature_id": simFactors.index,
        "biased_sim": biased_sim
    })
    
    # Write the bias information to the output file.
    bias_df.to_csv(bias_output_file, sep="\t", index=False)
    print(f"Bias file written to: {bias_output_file}")
    return simFactors, bias_df

def main():
    # Define file paths (update these with your actual locations)
    simulation_file = "./data/simulation-gene-readCounts.txt"
    bias_output_file = "./output/mapping-bias-genes-simulation.txt"
    
    # Run the simulation correction and write the bias output.
    simFactors, bias_df = simulation_correction(simulation_file, bias_output_file, analysis_functions)

if __name__ == "__main__":
    main()