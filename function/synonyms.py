import pandas as pd

def fbgn_to_data(fbgn_id):
    """
    Given an FBgn identifier, this function reads the FlyBase synonyms file 
    (FlyBase Synonyms file: fb_synonym_fb_2024_05.tsv) and returns all other data 
    associated with the provided FBgn identifier.
    
    The file includes:
      - primary_FBid: Primary FlyBase identifier for the object.
      - organism_abbreviation: Species abbreviation.
      - current_symbol: Current symbol used in FlyBase for the object.
      - current_fullname: Current full name used in FlyBase for the object.
      - fullname_synonym(s): Non-current full name(s) (pipe separated).
      - symbol_synonym(s): Non-current symbol(s) (pipe separated).
    
    Parameters:
        fbgn_id (str): The FBgn identifier to look up.
        
    Returns:
        dict: A dictionary containing all the data for the given FBgn identifier,
              excluding the 'primary_FBid' key. If no match is found, returns an empty dict.
    """
    file_path = "./database/fb_synonym_fb_2024_05.tsv"  # Update the path if necessary
    
    try:
        # Read the TSV file, skip the first 5 lines, and remove '#' from column names
        df = pd.read_csv(file_path, sep="\t", dtype=str, skiprows=5)
        df.columns = df.columns.str.replace('#', '')
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return {}
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return {}
    
    # Filter for rows where the primary_FBid matches the given fbgn_id
    match = df[df['primary_FBid'] == fbgn_id]
    
    if match.empty:
        return {}
    
    # There may be multiple rows, here we assume one record per FBgn.
    # Extract the first match and drop the primary_FBid column.
    record = match.iloc[0].to_dict()
    record.pop("primary_FBid", None)
    
    return record

def process_synonyms(data):
    """
    Process the synonym fields in the data dictionary, splitting the strings
    on the pipe character '|' and returning lists for each synonym type.
    
    Parameters:
        data (dict): The dictionary containing FlyBase data.
        
    Returns:
        dict: A new dictionary with the same keys as data, but with the 
              'fullname_synonym(s)' and 'symbol_synonym(s)' fields converted to lists.
    """
    processed_data = data.copy()
    
    if "fullname_synonym(s)" in processed_data and processed_data["fullname_synonym(s)"]:
        processed_data["fullname_synonym(s)"] = [syn.strip() for syn in processed_data["fullname_synonym(s)"].split("|")]
    
    if "symbol_synonym(s)" in processed_data and processed_data["symbol_synonym(s)"]:
        processed_data["symbol_synonym(s)"] = [syn.strip() for syn in processed_data["symbol_synonym(s)"].split("|")]
    
    return processed_data

def main():
    # Test the function with a sample FBgn identifier.
    test_fbgn = "FBgn0000303"  # Replace with an actual FBgn id for a real test
    data = fbgn_to_data(test_fbgn)
    
    if data:
        data = process_synonyms(data)
        print(f"Data for {test_fbgn}:")
        for key, value in data.items():
            # If the value is a list (i.e., synonyms), print each entry on a new line
            if isinstance(value, list):
                print(f"{key}:")
                for item in value:
                    print(f"  - {item}")
            else:
                print(f"{key}: {value}")
    else:
        print(f"No matching data found for FBgn id: {test_fbgn}")

if __name__ == "__main__":
    main()