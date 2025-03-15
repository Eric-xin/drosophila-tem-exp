import requests
import re

def fetch_sequence_location(species, location, strand='minus', padding=100):
    url = f"https://api.flybase.org/api/v1.0/sequence/region/{species}/{location}"
    params = {
        'strand': strand,
        'padding': padding
    }
    headers = {
        'accept': 'application/json'
    }
    
    response = requests.get(url, headers=headers, params=params)
    
    if response.status_code == 200:
        data = response.json()
        if "resultset" in data and "result" in data["resultset"]:
            return data["resultset"]["result"][0]["sequence"]
    else:
        response.raise_for_status()

# def convert_location_format(location_str):
#     parts = location_str.split('_')
#     return f"{parts[0]}:{parts[1]}..{parts[2]}"

def convert_location_format(location_str):
    parts = location_str.split('_')
    start = parts[1].split('.')[0]
    end = parts[2].split('.')[0]
    return f"{parts[0]}:{start}..{end}"

def is_location_sequence(seq_str):
    # Define the two patterns:
    # Pattern 1: Numeric chromosome with arm (e.g., "12L_123456_123789")
    pattern1 = re.compile(r"^\d+[LR]_\d+_\d+$")
    # Pattern 2: Alternate format (e.g., "X_7451422_7451891.0")
    pattern2 = re.compile(r"^[A-Z]+_\d+_\d+(?:\.0)?$")
    
    # Split the input by comma and remove extra spaces
    locations = [s.strip() for s in seq_str.split(',')]
    
    # Check each location string against the patterns
    for loc in locations:
        if not (pattern1.match(loc) or pattern2.match(loc)):
            return False
    return True

def v_location_format(location_str):
    # First, try the original pattern: e.g. "12L_123456_123789"
    match = re.match(r"(\d+)([LR])_(\d+)_(\d+)", location_str)
    if match:
        chrom, arm, start, end = match.groups()
        return f"{chrom}{arm}:{start}..{end}"
    
    # Next, try a pattern for locations like "X_7451422_7451891.0"
    match = re.match(r"([A-Z]+)_(\d+)_(\d+)(?:\.0)?$", location_str)
    if match:
        chrom, start, end = match.groups()
        return f"{chrom}:{start}..{end}"
    
    raise ValueError(f"Invalid location string: {location_str}")

# Example usage
if __name__ == "__main__":
    species = "dmel"
    # location = "2L:100000..101000"
    location = v_location_format("3R_6286767_6287794.0")
    data = fetch_sequence_location(species, location)
    print(data)