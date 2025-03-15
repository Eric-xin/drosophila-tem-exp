from flask import Flask, render_template, request, redirect, url_for, flash
import pandas as pd
import requests
import json
from werkzeug.utils import secure_filename
import os

app = Flask(__name__)
app.secret_key = 'your_secret_key'
app.config['UPLOAD_FOLDER'] = 'uploads'
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# --- API Functions ---
def fetch_sequence_id(gene_id):
    url = f"https://api.flybase.org/api/v1.0/sequence/id/{gene_id}"
    headers = {"accept": "application/json"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        try:
            data = response.json()
        except Exception:
            return "Error decoding JSON response."
        if "resultset" in data and "result" in data["resultset"]:
            return data["resultset"]["result"][0]["sequence"]
        else:
            return "No sequence found for the given gene ID."
    else:
        return f"Error: {response.status_code}"

def fetch_sequence_location(species, location, strand='minus', padding=100):
    species = species or "dmel"
    if ".." in location:
        location = location.replace("..", "-")
    url = f"https://api.flybase.org/api/v1.0/sequence/region/{species}/{location}"
    params = {'strand': strand, 'padding': padding}
    headers = {'accept': 'application/json'}
    response = requests.get(url, headers=headers, params=params)
    if response.status_code == 200:
        try:
            data = response.json()
        except Exception:
            return "Error decoding JSON response."
        if "resultset" in data and "result" in data["resultset"]:
            return data["resultset"]["result"][0]["sequence"]
        else:
            return "No sequence returned for this region."
    else:
        response.raise_for_status()

def fetch_id_abt(gene_id):
    url = f"https://api.flybase.org/api/v1.0/gene/summaries/auto/{gene_id}"
    headers = {"accept": "application/json"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        try:
            data = response.json()
        except Exception:
            return "Error decoding JSON response."
        if "resultset" in data and "result" in data["resultset"]:
            return data["resultset"]["result"][0]["summary"]
        else:
            return "No summary found for the given gene ID."
    else:
        return f"Error: {response.status_code}"

def format_sequence(seq):
    html = '<div style="white-space: normal;">'
    for char in seq:
        if char.upper() == "A":
            html += f'<span style="color: red; font-weight: bold;">{char}</span>'
        elif char.upper() == "T":
            html += f'<span style="color: blue; font-weight: bold;">{char}</span>'
        elif char.upper() == "C":
            html += f'<span style="color: green; font-weight: bold;">{char}</span>'
        elif char.upper() == "G":
            html += f'<span style="color: orange; font-weight: bold;">{char}</span>'
        else:
            html += char
    html += "</div>"
    return html

# Global variable to hold the DataFrame (for demo purposes)
df_global = None

# --- Preload TSV File ---
DEFAULT_TSV = os.path.join("data", "processed_allele_counts.tsv")
if os.path.exists(DEFAULT_TSV):
    try:
        df = pd.read_csv(DEFAULT_TSV, sep="\t")
        for col in ['start', 'end']:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
        for col in ['temp13_avg', 'temp18_avg', 'temp23_avg', 'temp29_avg']:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        df_global = df.copy()
        print("TSV file preloaded successfully.")
    except Exception as e:
        print("Error preloading TSV file:", e)
else:
    print(f"Default TSV file not found at: {DEFAULT_TSV}")

# --- Routes ---
@app.route('/')
def index():
    # If a TSV file was preloaded, show the table directly.
    if df_global is not None:
        table_html = df_global.to_html(classes='table table-striped', index=False)
        return render_template('data.html', table_html=table_html)
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload():
    if 'tsvfile' not in request.files:
        flash('No file part')
        return redirect(url_for('index'))
    file = request.files['tsvfile']
    if file.filename == '':
        flash('No selected file')
        return redirect(url_for('index'))
    filename = secure_filename(file.filename)
    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)
    try:
        df = pd.read_csv(filepath, sep="\t")
        for col in ['start', 'end']:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
        for col in ['temp13_avg', 'temp18_avg', 'temp23_avg', 'temp29_avg']:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        global df_global
        df_global = df.copy()  # Save globally for later access
        table_html = df.to_html(classes='table table-striped', index=False)
        return render_template('data.html', table_html=table_html)
    except Exception as e:
        flash(f"Error processing file: {e}")
        return redirect(url_for('index'))

@app.route('/detail/<int:row_id>')
def detail(row_id):
    global df_global
    if df_global is None or row_id < 0 or row_id >= len(df_global):
        flash("Invalid record selected.")
        return redirect(url_for('index'))
    row_data = df_global.iloc[row_id].to_dict()
    fbgn = str(row_data.get("fbgn", "")).strip()
    if fbgn and fbgn.startswith("FBgn"):
        row_data["Gene Summary"] = fetch_id_abt(fbgn)
        seq = fetch_sequence_id(fbgn)
        row_data["Gene Sequence"] = format_sequence(seq) if isinstance(seq, str) else seq
    return render_template('detail.html', row_data=row_data)

if __name__ == '__main__':
    app.run(debug=True)