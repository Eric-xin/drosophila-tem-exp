import pandas as pd

def fetch_flybase_record(fbgn_id):
    """
    Retrieve a FlyBase record for the given FBgn identifier from the synonyms file
    (fb_synonym_fb_2024_05.tsv) and remap its keys to concise names.

    Original keys and their new mappings:
      - primary_FBid: (removed)
      - organism_abbreviation -> org
      - current_symbol -> symbol
      - current_fullname -> name
      - fullname_synonym(s) -> name_synonyms
      - symbol_synonym(s) -> symbol_synonyms

    Parameters:
        fbgn_id (str): The FBgn identifier to look up.

    Returns:
        dict: A dictionary with concise keys and synonym fields converted to lists.
              If no match is found, returns an empty dict.
    """
    synonyms_filepath = "./database/fb_synonym_fb_2024_05.tsv"  # Update the path if necessary
    
    try:
        # Read the TSV file, skipping header lines, and clean column names
        df = pd.read_csv(synonyms_filepath, sep="\t", dtype=str, skiprows=5)
        df.columns = df.columns.str.replace('#', '')
    except FileNotFoundError:
        print(f"Error: File '{synonyms_filepath}' not found.")
        return {}
    except Exception as error:
        print(f"An error occurred while reading the file: {error}")
        return {}
    
    # Filter the dataframe for the matching FBgn id
    matching_records = df[df['primary_FBid'] == fbgn_id]
    if matching_records.empty:
        return {}
    
    # Extract the first matching record and convert to dictionary
    record = matching_records.iloc[0].to_dict()
    record.pop("primary_FBid", None)  # Remove the original primary key
    
    # Remap the remaining keys to simpler names
    key_mapping = {
        "organism_abbreviation": "org",
        "current_symbol": "symbol",
        "current_fullname": "name",
        "fullname_synonym(s)": "name_synonyms",
        "symbol_synonym(s)": "symbol_synonyms"
    }
    
    concise_record = {}
    for old_key, new_key in key_mapping.items():
        value = record.get(old_key, "")
        # Process synonym fields by splitting on '|'
        if new_key in ["name_synonyms", "symbol_synonyms"] and value:
            concise_record[new_key] = [item.strip() for item in value.split("|")]
        else:
            concise_record[new_key] = value
    
    return concise_record

def main():
    # Test the function with a sample FBgn identifier.
    sample_fbgn = "FBgn0000303"  # Replace with an actual FBgn id for a real test
    flybase_data = fetch_flybase_record(sample_fbgn)
    
    if flybase_data:
        print(f"Data for {sample_fbgn}:")
        for key, value in flybase_data.items():
            if isinstance(value, list):
                print(f"{key}:")
                for item in value:
                    print(f"  - {item}")
            else:
                print(f"{key}: {value}")
    else:
        print(f"No matching data found for FBgn id: {sample_fbgn}")

if __name__ == "__main__":
    main()