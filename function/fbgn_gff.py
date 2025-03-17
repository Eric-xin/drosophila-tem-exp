import gzip
import pandas as pd

def search_fbgn(gff_file, fbgn_id):
    """
    Search for and extract all GFF entries containing the given FBgn ID.
    Returns a pandas DataFrame with columns: seqid, source, feature, start, end,
    score, strand, phase, and attributes.
    """
    # Determine whether to use gzip.open or regular open based on file extension
    open_func = gzip.open if gff_file.endswith('.gz') else open
    
    data = []
    with open_func(gff_file, 'rt') as f:
        for line in f:
            # Skip comment lines
            if line.startswith("#"):
                continue
            if fbgn_id in line:
                # Split line into GFF columns (tab-separated)
                cols = line.strip().split("\t")
                if len(cols) >= 9:
                    data.append(cols[:9])
                    
    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    df = pd.DataFrame(data, columns=columns)
    return df

def main():
    # Define test input parameters
    gff_file = "./database/dmel-all-r6.62.gff.gz"  # Replace with your actual file path
    fbgn_id = "FBgn0034730"             # Replace with your FBgn ID of interest
    output_file = "result.csv"          # Define the output CSV file path
    
    # Run the search function and get a DataFrame
    df = search_fbgn(gff_file, fbgn_id)
    
    if not df.empty:
        # Save the DataFrame to the output CSV file
        df.to_csv(output_file, index=False)
        print(f"Found {len(df)} matching entries. Results written to {output_file}.")
    else:
        print(f"FBgn ID {fbgn_id} not found in {gff_file}.")

if __name__ == '__main__':
    main()