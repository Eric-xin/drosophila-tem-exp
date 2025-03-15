import pandas as pd

def cg_to_fbgn(cg_id):
    """
    Given a CG identifier, this function reads the fixed FlyBase annotation file
    (./database/fbgn_annotation_ID.tsv.gz) and returns the corresponding primary
    and secondary FBgn identifiers.
    
    Parameters:
        cg_id (str): The CG identifier to look up.
        
    Returns:
        tuple: A tuple (primary_fbgn, secondary_fbgn) where:
            - primary_fbgn is a list of primary FBgn identifiers.
            - secondary_fbgn is a list of secondary FBgn identifiers.
        Both lists will be empty if no match is found.
    """
    file_path = "./database/fbgn_annotation_ID_fb_2024_05.tsv"
    # Read the compressed TSV file. The 'compression' parameter is set to infer automatically.
    # df = pd.read_csv(file_path, sep="\t", compression='infer', dtype=str)
    df = pd.read_csv(file_path, sep="\t", dtype=str)
    
    def in_secondary_ids(cell):
        if pd.isna(cell):
            return False
        return cg_id in [x.strip() for x in cell.split(',')]
    
    # Create a mask checking if cg_id matches the primary annotation_ID or is present in the secondary annotation_ID(s)
    mask = (df['annotation_ID'] == cg_id) | (df['secondary_annotation_ID(s)'].apply(in_secondary_ids))
    
    if not mask.any():
        return ([], [])
    
    primary_fbgn = set()
    secondary_fbgn = set()
    
    for _, row in df[mask].iterrows():
        # Add the primary FBgn identifier.
        primary_fbgn.add(row['primary_FBgn#'])
        # If secondary FBgn identifiers are present, add them.
        if pd.notna(row['secondary_FBgn#(s)']):
            ids = [s.strip() for s in row['secondary_FBgn#(s)'].split(',')]
            secondary_fbgn.update(ids)
    
    return (list(primary_fbgn), list(secondary_fbgn))

def main():
    # Test the function with a sample CG identifier.
    test_cg = "CG12345"  # Replace with an actual CG id for a real test.
    primary, secondary = cg_to_fbgn(test_cg)
    
    if primary or secondary:
        print(f"CG id: {test_cg}")
        print("Primary FBgn id(s):", ", ".join(primary) if primary else "None")
        print("Secondary FBgn id(s):", ", ".join(secondary) if secondary else "None")
    else:
        print(f"No matching FBgn found for CG id: {test_cg}")

if __name__ == "__main__":
    main()