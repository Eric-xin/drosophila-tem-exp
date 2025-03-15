#!/usr/bin/env python3
import pandas as pd
import re
from tqdm import tqdm

#########################################
# Helper functions for FlyBase lookups  #
#########################################

def v_location_format(location_str):
    """
    Reformat a location string.
    
    First, try the original pattern: e.g. "12L_123456_123789"
    Next, try a pattern for locations like "X_7451422_7451891.0"
    """
    match = re.match(r"(\d+)([LR])_(\d+)_(\d+)", location_str)
    if match:
        chrom, arm, start, end = match.groups()
        return f"{chrom}{arm}:{start}..{end}"
    
    match = re.match(r"([A-Z]+)_(\d+)_(\d+)(?:\.0)?$", location_str)
    if match:
        chrom, start, end = match.groups()
        return f"{chrom}:{start}..{end}"
    
    raise ValueError(f"Invalid location string: {location_str}")

def parse_location(formatted_location):
    """
    Given a formatted location (e.g. "12L:123456..123789"),
    return the arm, start, and end.
    """
    try:
        arm_part, coords = formatted_location.split(":")
        start, end = coords.split("..")
        return arm_part, start, end
    except Exception:
        return "", "", ""

def fetch_flybase_record(fbgn_id):
    """
    Retrieve a FlyBase record for the given FBgn identifier
    and remap its keys to concise names.
    
    Mapped keys:
      - organism_abbreviation -> org
      - current_symbol        -> symbol
      - current_fullname      -> name
      - fullname_synonym(s)   -> name_synonyms
      - symbol_synonym(s)     -> symbol_synonyms
    """
    synonyms_filepath = "./database/fb_synonym_fb_2024_05.tsv"  # Update path if needed
    try:
        df = pd.read_csv(synonyms_filepath, sep="\t", dtype=str, skiprows=5)
        df.columns = df.columns.str.replace('#', '')
    except FileNotFoundError:
        print(f"Error: File '{synonyms_filepath}' not found.")
        return {}
    except Exception as error:
        print(f"Error reading FlyBase synonyms file: {error}")
        return {}
    
    matching_records = df[df['primary_FBid'] == fbgn_id]
    if matching_records.empty:
        return {}
    
    record = matching_records.iloc[0].to_dict()
    record.pop("primary_FBid", None)
    
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
        if new_key in ["name_synonyms", "symbol_synonyms"]:
            if isinstance(value, str) and value:
                concise_record[new_key] = [item.strip() for item in value.split("|")]
            else:
                concise_record[new_key] = []
        else:
            concise_record[new_key] = value
    return concise_record

def fbgn_to_reference(fbgn_id):
    """
    Retrieve reference location data for a given FBgn identifier.
    Expected columns in the TSV file:
      - FBgn_ID, Arm/Scaffold, Location, Strand
    Splits the Location field (e.g. "18705501..18732344") into Start and End.
    """
    file_path = "./database/dmel_paralogs_fb_2024_05.tsv"  # Update path if needed
    try:
        df = pd.read_csv(file_path, sep="\t", dtype=str)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return {}
    except Exception as e:
        print(f"Error reading reference file: {e}")
        return {}
    
    record = df[df['FBgn_ID'] == fbgn_id]
    if record.empty:
        return {}
    
    row = record.iloc[0]
    arm = row.get('Arm/Scaffold', "")
    location = row.get('Location', "")
    strand = row.get('Strand', "")
    
    try:
        start, end = location.split("..")
    except Exception:
        start, end = "", ""
    
    return {"Arm": arm, "Start": start, "End": end, "Strand": strand}

def cg_to_fbgn(cg_id):
    """
    Convert a CG identifier to its corresponding FBgn identifiers.
    
    Returns a tuple (primary_fbgn, secondary_fbgn) where each is a list.
    """
    file_path = "./database/fbgn_annotation_ID_fb_2024_05.tsv"
    try:
        df = pd.read_csv(file_path, sep="\t", dtype=str)
    except Exception as e:
        print(f"Error reading CG-to-FBgn file: {e}")
        return ([], [])
    
    def in_secondary_ids(cell):
        if pd.isna(cell):
            return False
        return cg_id in [x.strip() for x in cell.split(',')]
    
    mask = (df['annotation_ID'] == cg_id) | (df['secondary_annotation_ID(s)'].apply(in_secondary_ids))
    if not mask.any():
        return ([], [])
    
    primary_fbgn = set()
    secondary_fbgn = set()
    for _, row in df[mask].iterrows():
        primary_fbgn.add(row['primary_FBgn#'])
        if pd.notna(row.get('secondary_FBgn#(s)', "")):
            ids = [s.strip() for s in row['secondary_FBgn#(s)'].split(',')]
            secondary_fbgn.update(ids)
    
    return (list(primary_fbgn), list(secondary_fbgn))

#########################################
# Main processing script                #
#########################################

def main():
    # Read the allele counts table (adjust path and separator as needed)
    input_file = "./data/gene-readCounts-F1.txt"
    try:
        df = pd.read_csv(input_file, sep="\t", dtype=str)
    except Exception as e:
        print(f"Error reading allele counts table: {e}")
        return
    
    # Create new columns for annotations and location info.
    df["fbgn"] = ""
    df["current_name"] = ""
    df["arm"] = ""
    df["start"] = ""
    df["end"] = ""
    df["strand"] = ""
    
    valid_indices = []
    invalid_indices = []
    
    # Process each row by examining the feature_id.
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing rows"):
        feature_id = row["feature_id"]
        fbgn = ""
        # If feature_id is an FBgn identifier, process accordingly.
        if feature_id.startswith("FBgn"):
            fbgn = feature_id
        # If feature_id is a CG identifier, convert to FBgn.
        elif feature_id.startswith("CG"):
            primary, secondary = cg_to_fbgn(feature_id)
            if primary:
                fbgn = primary[0]
            elif secondary:
                fbgn = secondary[0]
        # If neither FBgn nor CG, try to treat feature_id as a location string.
        else:
            try:
                formatted = v_location_format(feature_id)
                arm_parsed, start_parsed, end_parsed = parse_location(formatted)
                df.at[idx, "arm"] = arm_parsed
                df.at[idx, "start"] = start_parsed
                df.at[idx, "end"] = end_parsed
                valid_indices.append(idx)
                continue  # Skip further processing for this row.
            except ValueError:
                invalid_indices.append(idx)
                continue
        
        # If a valid FBgn was determined, annotate the row.
        if fbgn:
            df.at[idx, "fbgn"] = fbgn
            record = fetch_flybase_record(fbgn)
            if record:
                df.at[idx, "current_name"] = record.get("name", "")
            loc = fbgn_to_reference(fbgn)
            if loc:
                try:
                    loc_string = f"{loc['Arm']}_{loc['Start']}_{loc['End']}"
                    formatted = v_location_format(loc_string)
                    arm_parsed, start_parsed, end_parsed = parse_location(formatted)
                    df.at[idx, "arm"] = arm_parsed
                    df.at[idx, "start"] = start_parsed
                    df.at[idx, "end"] = end_parsed
                except Exception as e:
                    df.at[idx, "arm"] = loc.get("Arm", "")
                    df.at[idx, "start"] = loc.get("Start", "")
                    df.at[idx, "end"] = loc.get("End", "")
                df.at[idx, "strand"] = loc.get("Strand", "")
            valid_indices.append(idx)
        else:
            invalid_indices.append(idx)
    
    ##########################################################
    # Calculate overall expression averages per temperature  #
    ##########################################################
    # Temperature groups defined by replicate numbers:
    # 1-3 -> 13°C, 4-6 -> 18°C, 7-9 -> 23°C, 10-12 -> 29°C.
    temp_groups = {
        "temp13": [1, 2, 3],
        "temp18": [4, 5, 6],
        "temp23": [7, 8, 9],
        "temp29": [10, 11, 12]
    }
    
    # Get all replicate columns (they start with "rep")
    rep_cols = [col for col in df.columns if col.startswith("rep")]
    
    # For each temperature group, select all replicate columns whose replicate number is in the group
    for temp, rep_nums in temp_groups.items():
        temp_cols = []
        for col in rep_cols:
            m = re.search(r"rep(\d+)_", col)
            if m:
                rep_num = int(m.group(1))
                if rep_num in rep_nums:
                    temp_cols.append(col)
        # Calculate the overall average (regardless of allele) for this temperature group
        df[temp + "_avg"] = df[temp_cols].apply(lambda row: pd.to_numeric(row, errors='coerce').mean(), axis=1)
    
    #########################################
    # Remove original columns               #
    #########################################
    # Drop the original feature_id and all replicate columns.
    cols_to_drop = [col for col in df.columns if col == "feature_id" or col.startswith("rep")]
    
    # For the processed (valid) rows:
    df_valid = df.loc[valid_indices].copy()
    df_valid.drop(columns=cols_to_drop, inplace=True)
    
    # For unfiltered rows (where conversion failed), keep them as-is.
    df_invalid = df.loc[invalid_indices].copy()
    
    #########################################
    # Output results                        #
    #########################################
    df_valid.to_csv("./output/processed_allele_counts.tsv", sep="\t", index=False)
    df_invalid.to_csv("./output/unfiltered_allele_counts.tsv", sep="\t", index=False)
    
    print("Processing complete.")
    print("  - Processed data saved to './output/processed_allele_counts.tsv'")
    print("  - Unfiltered rows saved to './output/unfiltered_allele_counts.tsv'")

if __name__ == "__main__":
    main()