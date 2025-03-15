import pandas as pd

def fbgn_to_reference(fbgn_id):
    """
    Retrieve reference data for a given FlyBase FBgn identifier from a TSV file.
    
    The TSV file is expected to have at least the following columns:
      - FBgn_ID: Current FlyBase identifier (FBgn#) of the D. melanogaster gene.
      - GeneSymbol: Current FlyBase gene symbol.
      - Arm/Scaffold: Arm or scaffold where the gene is localized.
      - Location: Genomic coordinates of the gene on the arm, e.g., "18705501..18732344".
      - Strand: Orientation of the gene ('1' for positive, '-1' for negative).
    
    Parameters:
        fbgn_id (str): The FBgn identifier to look up.
        
    Returns:
        dict: A dictionary with keys "Arm/Scaffold", "Start", "End", and "Strand" for the given FBgn identifier.
              If no match is found, returns an empty dict.
    """
    file_path = "./database/dmel_paralogs_fb_2024_05.tsv"  # Update the path as needed
    
    try:
        # Read the TSV file (assuming the first row contains headers)
        df = pd.read_csv(file_path, sep="\t", dtype=str)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return {}
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return {}
    
    # Filter the DataFrame for the matching FBgn_ID
    record = df[df['FBgn_ID'] == fbgn_id]
    
    if record.empty:
        return {}
    
    # Extract the first matching record and the relevant columns
    row = record.iloc[0]
    arm = row['Arm/Scaffold']
    location = row['Location']
    strand = row['Strand']
    
    # Split the location field on '..' into start and end positions
    try:
        start, end = location.split("..")
    except Exception as e:
        print(f"Error parsing location '{location}': {e}")
        start, end = None, None
    
    return {
        "Arm": arm,
        "Start": start,
        "End": end,
        "Strand": strand
    }

def main():
    # Example FBgn identifier; replace with an actual FBgn ID as needed.
    test_fbgn = "FBgn0000303"
    data = fbgn_to_reference(test_fbgn)
    
    if data:
        print(f"Reference data for {test_fbgn}:")
        for key, value in data.items():
            print(f"{key}: {value}")
    else:
        print(f"No reference data found for FBgn ID: {test_fbgn}")

if __name__ == "__main__":
    main()