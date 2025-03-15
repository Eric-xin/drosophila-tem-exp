import requests

def fetch_id_abt(gene_id):
    url = f"https://api.flybase.org/api/v1.0/gene/summaries/auto/{gene_id}"
    headers = {"accept": "application/json"}
    
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        data = response.json()
        if "resultset" in data and "result" in data["resultset"]:
            return data["resultset"]["result"][0]["summary"]
        else:
            return "No sequence found for the given gene ID."
    else:
        return f"Error: {response.status_code}"

# Example usage
if __name__ == "__main__":
    gene_id = "FBgn0027621"  # Replace with the desired gene ID
    sequence = fetch_id_abt(gene_id)
    print(sequence)