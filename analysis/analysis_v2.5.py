#!/usr/bin/env python3
import pandas as pd
import re
import concurrent.futures
from tqdm import tqdm

#########################################
# Global variables for loaded databases #
#########################################
DF_SYNONYMS = None
DF_REFERENCE = None
DF_CG = None

#########################################
# Worker initializer                   #
#########################################
def init_worker(synonyms_df, reference_df, cg_df):
    global DF_SYNONYMS, DF_REFERENCE, DF_CG
    DF_SYNONYMS = synonyms_df
    DF_REFERENCE = reference_df
    DF_CG = cg_df

#########################################
# Helper functions for processing rows   #
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
    from the preloaded DF_SYNONYMS and remap its keys.
    """
    global DF_SYNONYMS
    if DF_SYNONYMS is None:
        return {}
    matching = DF_SYNONYMS[DF_SYNONYMS['primary_FBid'] == fbgn_id]
    if matching.empty:
        return {}
    record = matching.iloc[0].to_dict()
    record.pop("primary_FBid", None)
    key_mapping = {
        "organism_abbreviation": "org",
        "current_symbol": "symbol",
        "current_fullname": "name",
        "fullname_synonym(s)": "name_synonyms",
        "symbol_synonym(s)": "symbol_synonyms"
    }
    concise = {}
    for old_key, new_key in key_mapping.items():
        value = record.get(old_key, "")
        if new_key in ["name_synonyms", "symbol_synonyms"]:
            if isinstance(value, str) and value:
                concise[new_key] = [item.strip() for item in value.split("|")]
            else:
                concise[new_key] = []
        else:
            concise[new_key] = value
    return concise

def fbgn_to_reference(fbgn_id):
    """
    Retrieve reference location data for a given FBgn identifier
    from the preloaded DF_REFERENCE.
    """
    global DF_REFERENCE
    if DF_REFERENCE is None:
        return {}
    rec = DF_REFERENCE[DF_REFERENCE['FBgn_ID'] == fbgn_id]
    if rec.empty:
        return {}
    row = rec.iloc[0]
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
    Convert a CG identifier to FBgn identifiers using the preloaded DF_CG.
    Returns a tuple (primary_fbgn, secondary_fbgn) where each is a list.
    """
    global DF_CG
    if DF_CG is None:
        return ([], [])
    def in_secondary(cell):
        if pd.isna(cell):
            return False
        return cg_id in [x.strip() for x in cell.split(',')]
    mask = (DF_CG['annotation_ID'] == cg_id) | (DF_CG['secondary_annotation_ID(s)'].apply(in_secondary))
    if not mask.any():
        return ([], [])
    primary = set()
    secondary = set()
    for _, row in DF_CG[mask].iterrows():
        primary.add(row['primary_FBgn#'])
        if pd.notna(row.get('secondary_FBgn#(s)', "")):
            ids = [s.strip() for s in row['secondary_FBgn#(s)'].split(',')]
            secondary.update(ids)
    return (list(primary), list(secondary))

def process_row(item):
    """
    Process a single row.
    item: tuple (index, row_dict)
    Returns a tuple (index, updated_row_dict, valid_flag)
    """
    idx, row = item
    new_row = row.copy()
    fbgn = ""
    valid = False
    feature_id = row.get("feature_id", "")
    # If FBgn:
    if feature_id.startswith("FBgn"):
        fbgn = feature_id
    # If CG:
    elif feature_id.startswith("CG"):
        primary, secondary = cg_to_fbgn(feature_id)
        if primary:
            fbgn = primary[0]
        elif secondary:
            fbgn = secondary[0]
    # Otherwise, try to treat as a location string.
    else:
        try:
            formatted = v_location_format(feature_id)
            arm_parsed, start_parsed, end_parsed = parse_location(formatted)
            new_row["arm"] = arm_parsed
            new_row["start"] = start_parsed
            new_row["end"] = end_parsed
            valid = True
            return (idx, new_row, valid)
        except ValueError:
            return (idx, new_row, False)
    # If we got a valid FBgn, update annotations.
    if fbgn:
        new_row["fbgn"] = fbgn
        rec = fetch_flybase_record(fbgn)
        if rec:
            new_row["current_name"] = rec.get("name", "")
        loc = fbgn_to_reference(fbgn)
        if loc:
            try:
                loc_string = f"{loc['Arm']}_{loc['Start']}_{loc['End']}"
                formatted = v_location_format(loc_string)
                arm_parsed, start_parsed, end_parsed = parse_location(formatted)
                new_row["arm"] = arm_parsed
                new_row["start"] = start_parsed
                new_row["end"] = end_parsed
            except Exception:
                new_row["arm"] = loc.get("Arm", "")
                new_row["start"] = loc.get("Start", "")
                new_row["end"] = loc.get("End", "")
            new_row["strand"] = loc.get("Strand", "")
        valid = True
    return (idx, new_row, valid)

#########################################
# Main processing script                #
#########################################

def main():
    # Load the allele counts table.
    input_file = "./data/gene-readCounts-F1.txt"
    try:
        df = pd.read_csv(input_file, sep="\t", dtype=str)
    except Exception as e:
        print(f"Error reading allele counts table: {e}")
        return

    # Preload database files once.
    try:
        df_syn = pd.read_csv("./database/fb_synonym_fb_2024_05.tsv", sep="\t", dtype=str, skiprows=5)
        df_syn.columns = df_syn.columns.str.replace('#', '')
    except Exception as e:
        print(f"Error loading FlyBase synonyms file: {e}")
        return
    try:
        df_ref = pd.read_csv("./database/dmel_paralogs_fb_2024_05.tsv", sep="\t", dtype=str)
    except Exception as e:
        print(f"Error loading reference file: {e}")
        return
    try:
        df_cg = pd.read_csv("./database/fbgn_annotation_ID_fb_2024_05.tsv", sep="\t", dtype=str)
    except Exception as e:
        print(f"Error loading CG-to-FBgn file: {e}")
        return

    # Prepare lists to hold indices for valid and invalid rows.
    valid_indices = []
    invalid_indices = []

    # Use a ProcessPoolExecutor to process rows in parallel.
    items = list(df.iterrows())
    processed_results = []
    with concurrent.futures.ProcessPoolExecutor(
            initializer=init_worker,
            initargs=(df_syn, df_ref, df_cg)
        ) as executor:
        # Use tqdm for progress.
        for result in tqdm(executor.map(process_row, items), total=len(items), desc="Processing rows"):
            processed_results.append(result)

    # Reassemble results.
    new_rows = {}
    for idx, new_row, valid in processed_results:
        new_rows[idx] = new_row
        if valid:
            valid_indices.append(idx)
        else:
            invalid_indices.append(idx)
    # Replace rows in the original DataFrame with processed ones.
    df_processed = pd.DataFrame.from_dict(new_rows, orient="index")
    df_processed = df_processed.sort_index()

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
    # All replicate columns (assumed to start with "rep")
    rep_cols = [col for col in df_processed.columns if col.startswith("rep")]
    for temp, rep_nums in temp_groups.items():
        temp_cols = []
        for col in rep_cols:
            m = re.search(r"rep(\d+)_", col)
            if m:
                rep_num = int(m.group(1))
                if rep_num in rep_nums:
                    temp_cols.append(col)
        df_processed[temp + "_avg"] = df_processed[temp_cols].apply(
            lambda row: pd.to_numeric(row, errors='coerce').mean(), axis=1
        )
    
    #########################################
    # Remove original columns               #
    #########################################
    # Drop the original feature_id and replicate columns.
    cols_to_drop = [col for col in df_processed.columns if col == "feature_id" or col.startswith("rep")]
    df_processed.drop(columns=cols_to_drop, inplace=True)

    # Split into valid and invalid sets.
    df_valid = df_processed.loc[valid_indices].copy()
    df_invalid = df_processed.loc[invalid_indices].copy()

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