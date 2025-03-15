import sys
import os
import json
import pandas as pd
import requests
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
    QFileDialog, QComboBox, QPushButton, QHBoxLayout, QLabel, QDialog, QScrollArea,
    QGroupBox, QMessageBox, QHeaderView, QGridLayout, QTextEdit
)
from PyQt5.QtGui import QFont, QColor
from PyQt5.QtCore import Qt

# --- FlyBase API functions with improved error handling ---

def fetch_sequence_id(gene_id):
    """Fetch gene sequence for a given FBgn gene_id."""
    url = f"https://api.flybase.org/api/v1.0/sequence/id/{gene_id}"
    headers = {"accept": "application/json"}
    
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        try:
            data = response.json()
        except Exception as e:
            return "Error decoding JSON response."
        if "resultset" in data and "result" in data["resultset"]:
            return data["resultset"]["result"][0]["sequence"]
        else:
            return "No sequence found for the given gene ID."
    else:
        return f"Error: {response.status_code}"

def fetch_sequence_location(species, location, strand='minus', padding=100):
    """Fetch sequence by location if no FBgn ID is present.
       The location is expected in a format like '2L:10710118-10710729'."""
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
        try:
            data = response.json()
        except Exception as e:
            return "Error decoding JSON response."
        if "resultset" in data and "result" in data["resultset"]:
            return data["resultset"]["result"][0]["sequence"]
        else:
            return "No sequence returned for this region."
    else:
        response.raise_for_status()

def fetch_id_abt(gene_id):
    """Fetch gene summary (about) information for a given FBgn gene_id."""
    url = f"https://api.flybase.org/api/v1.0/gene/summaries/auto/{gene_id}"
    headers = {"accept": "application/json"}
    
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        try:
            data = response.json()
        except Exception as e:
            return "Error decoding JSON response."
        if "resultset" in data and "result" in data["resultset"]:
            return data["resultset"]["result"][0]["summary"]
        else:
            return "No summary found for the given gene ID."
    else:
        return f"Error: {response.status_code}"

# --- Utility function for formatting sequence with colors ---

def format_sequence(seq):
    """Return an HTML formatted sequence with stylized ATCG colors and wrapping."""
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

# --- Detail Dialog using QTextEdit for multi-line display ---

class DetailDialog(QDialog):
    def __init__(self, row_data, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Gene Detail")
        self.resize(700, 500)
        self.details_dict = dict(row_data)  # copy initial data

        main_layout = QVBoxLayout(self)
        
        # "Copy as JSON" button
        copy_button = QPushButton("Copy as JSON")
        copy_button.clicked.connect(self.copy_as_json)
        main_layout.addWidget(copy_button)
        
        # Scroll area for details
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        main_layout.addWidget(scroll)
        
        content = QWidget()
        scroll.setWidget(content)
        # Use a grid layout (2 columns) for the detail boxes
        grid = QGridLayout(content)
        grid.setSpacing(10)
        self.detail_layout = grid
        
        # Add initial data fields into the grid
        i = 0
        for key, value in row_data.items():
            widget = self.create_field_box(key, value)
            grid.addWidget(widget, i // 2, i % 2)
            i += 1
        
        # Fetch additional data
        fbgn = row_data.get("fbgn", "").strip()
        arm = row_data.get("arm", "").strip()
        start = row_data.get("start", "")
        end = row_data.get("end", "")
        
        if fbgn:
            # If FBgn ID is present, fetch both summary and sequence.
            summary = fetch_id_abt(fbgn)
            self.details_dict["Gene Summary"] = summary
            widget = self.create_field_box("Gene Summary", summary)
            grid.addWidget(widget, i // 2, i % 2)
            i += 1
            
            seq = fetch_sequence_id(fbgn)
            self.details_dict["Gene Sequence"] = seq
            widget = self.create_sequence_box("Gene Sequence", seq)
            grid.addWidget(widget, i // 2, i % 2)
            i += 1
        elif arm and start and end:
            # If no FBgn but location is available, only fetch sequence using location.
            location = f"{arm}:{start}-{end}"
            seq = fetch_sequence_location("dmel", location, strand=row_data.get("strand", "minus"))
            self.details_dict["Gene Sequence"] = seq
            widget = self.create_sequence_box(f"Gene Sequence (from {location})", seq)
            grid.addWidget(widget, i // 2, i % 2)
            i += 1
        # Otherwise, if no FBgn and no location, do not scrape additional information.
    
    def create_field_box(self, title, content):
        """Create a styled box for a single field that shows all text with wrapping using QTextEdit."""
        box = QGroupBox(title)
        layout = QVBoxLayout(box)
        text_edit = QTextEdit()
        text_edit.setPlainText(str(content))
        text_edit.setReadOnly(True)
        # Remove border for a cleaner look
        text_edit.setFrameStyle(0)
        # Ensure word wrap is enabled
        text_edit.setLineWrapMode(QTextEdit.WidgetWidth)
        layout.addWidget(text_edit)
        return box
    
    def create_sequence_box(self, title, sequence):
        """Create a styled box for gene sequence with ATCG styling and multi-line wrapping using QTextEdit."""
        box = QGroupBox(title)
        layout = QVBoxLayout(box)
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setFrameStyle(0)
        text_edit.setLineWrapMode(QTextEdit.WidgetWidth)
        if sequence and isinstance(sequence, str):
            text_edit.setHtml(format_sequence(sequence))
        else:
            text_edit.setPlainText("No sequence available.")
        layout.addWidget(text_edit)
        return box
    
    def copy_as_json(self):
        """Copy the details dictionary as JSON to the clipboard."""
        json_str = json.dumps(self.details_dict, indent=4)
        clipboard = QApplication.clipboard()
        clipboard.setText(json_str)
        QMessageBox.information(self, "Copied", "Data copied as JSON to clipboard.")

# --- Main Data Explorer Window ---
# (The rest of the code remains unchanged; see previous versions for file loading, sorting, and table display)

class DataExplorer(QMainWindow):
    def __init__(self, preload_file=None):
        super().__init__()
        self.setWindowTitle("Data Explorer")
        self.resize(1000, 600)
        self.df = None  # will hold the pandas DataFrame
        
        # Main widget and layout
        main_widget = QWidget(self)
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout(main_widget)
        
        # --- Top control panel ---
        control_layout = QHBoxLayout()
        
        # Button to load TSV file
        self.loadButton = QPushButton("Load TSV Data")
        self.loadButton.clicked.connect(self.open_file_dialog)
        control_layout.addWidget(self.loadButton)
        
        # Label for sort mode
        control_layout.addWidget(QLabel("Sort Mode:"))
        
        # ComboBox for selecting sort mode
        self.sortModeCombo = QComboBox()
        self.sortModeCombo.addItems([
            "Sort by Location",
            "Sort by Temperature Reading",
            "Sort by Temperature Ratio"
        ])
        self.sortModeCombo.currentIndexChanged.connect(self.on_sort_mode_change)
        control_layout.addWidget(self.sortModeCombo)
        
        # ComboBox for Temperature Reading sort mode
        self.tempSortCombo = QComboBox()
        self.tempSortCombo.addItems(["temp13_avg", "temp18_avg", "temp23_avg", "temp29_avg"])
        self.tempSortCombo.setEnabled(False)
        control_layout.addWidget(QLabel("Temperature Reading:"))
        control_layout.addWidget(self.tempSortCombo)
        
        # For Temperature Ratio, we need two different temperature selections:
        self.ratioNumCombo = QComboBox()
        self.ratioNumCombo.addItems(["temp13_avg", "temp18_avg", "temp23_avg", "temp29_avg"])
        self.ratioNumCombo.setEnabled(False)
        control_layout.addWidget(QLabel("Temperature Numerator:"))
        control_layout.addWidget(self.ratioNumCombo)
        
        self.ratioDenomCombo = QComboBox()
        self.ratioDenomCombo.addItems(["temp13_avg", "temp18_avg", "temp23_avg", "temp29_avg"])
        self.ratioDenomCombo.setEnabled(False)
        control_layout.addWidget(QLabel("Temperature Denom:"))
        control_layout.addWidget(self.ratioDenomCombo)
        
        # Button to apply sorting
        self.sortButton = QPushButton("Sort Data")
        self.sortButton.clicked.connect(self.sort_data)
        control_layout.addWidget(self.sortButton)
        
        main_layout.addLayout(control_layout)
        
        # --- Table Widget to show data ---
        self.table = QTableWidget(self)
        self.table.setColumnCount(0)
        self.table.setRowCount(0)
        self.table.setSortingEnabled(False)
        self.table.cellDoubleClicked.connect(self.table_double_click)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.setStyleSheet("""
            QTableWidget {
                background-color: #f9f9f9;
                font-size: 12pt;
            }
            QHeaderView::section {
                background-color: #0078d7;
                color: white;
                padding: 4px;
                font-size: 12pt;
            }
            QTableWidget::item:selected {
                background-color: #cce6ff;
            }
        """)
        main_layout.addWidget(self.table)
        
        # If a preload file was provided, load it immediately.
        if preload_file is not None:
            self.load_file(preload_file)
    
    def on_sort_mode_change(self, index):
        """Enable or disable the temperature combo boxes based on sort mode."""
        mode = self.sortModeCombo.currentText()
        if mode == "Sort by Temperature Ratio":
            self.ratioNumCombo.setEnabled(True)
            self.ratioDenomCombo.setEnabled(True)
            self.tempSortCombo.setEnabled(False)
        elif mode == "Sort by Temperature Reading":
            self.tempSortCombo.setEnabled(True)
            self.ratioNumCombo.setEnabled(False)
            self.ratioDenomCombo.setEnabled(False)
        else:
            self.tempSortCombo.setEnabled(False)
            self.ratioNumCombo.setEnabled(False)
            self.ratioDenomCombo.setEnabled(False)
    
    def open_file_dialog(self):
        """Open a file dialog and load the selected TSV file."""
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open TSV File", "", "TSV Files (*.tsv *.txt);;All Files (*)", options=options
        )
        if filename:
            self.load_file(filename)
    
    def load_file(self, filename):
        """Load the TSV file at the given filename into a DataFrame."""
        try:
            self.df = pd.read_csv(filename, sep="\t")
            # Ensure numeric columns are read as numbers
            for col in ['start', 'end', 'temp13_avg', 'temp18_avg', 'temp23_avg', 'temp29_avg']:
                if col in self.df.columns:
                    self.df[col] = pd.to_numeric(self.df[col], errors='coerce')
            self.populate_table(self.df)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading file: {e}")
    
    def populate_table(self, df):
        """Populate the QTableWidget with data from the DataFrame."""
        columns = list(df.columns)
        if "ratio" in df.columns and self.sortModeCombo.currentText() == "Sort by Temperature Ratio":
            if "ratio" not in columns:
                columns.append("ratio")
        
        self.table.setColumnCount(len(columns))
        self.table.setHorizontalHeaderLabels(columns)
        self.table.setRowCount(len(df))
        
        for row in range(len(df)):
            for col, col_name in enumerate(columns):
                value = df.iloc[row][col_name] if col_name in df.columns else ""
                if pd.isna(value):
                    cell_text = "N/A"
                else:
                    if isinstance(value, float):
                        cell_text = f"{value:.4f}"
                    else:
                        cell_text = str(value)
                item = QTableWidgetItem(cell_text)
                if isinstance(value, (int, float)) and not pd.isna(value):
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                else:
                    item.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)
                if pd.isna(value):
                    font = item.font()
                    font.setStrikeOut(True)
                    item.setFont(font)
                    item.setForeground(QColor('gray'))
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()
    
    def sort_data(self):
        """Sort the data based on the selected mode."""
        if self.df is None:
            return
        
        mode = self.sortModeCombo.currentText()
        df_sorted = self.df.copy()
        if mode == "Sort by Location":
            sort_cols = []
            if "arm" in df_sorted.columns:
                sort_cols.append("arm")
            if "start" in df_sorted.columns:
                sort_cols.append("start")
            if "end" in df_sorted.columns:
                sort_cols.append("end")
            df_sorted = df_sorted.sort_values(by=sort_cols, ascending=True)
        elif mode == "Sort by Temperature Reading":
            sort_col = self.tempSortCombo.currentText()
            df_sorted = df_sorted.sort_values(by=sort_col, ascending=False, na_position='last')
        elif mode == "Sort by Temperature Ratio":
            num_col = self.ratioNumCombo.currentText()
            denom_col = self.ratioDenomCombo.currentText()
            if num_col == denom_col:
                QMessageBox.warning(self, "Invalid Selection", "Please select two different temperature columns for the ratio.")
                return
            def compute_ratio(row):
                try:
                    numerator = float(row.get(num_col, 0))
                    denominator = float(row.get(denom_col, 0))
                    return numerator / denominator if denominator != 0 else float('nan')
                except:
                    return float('nan')
            df_sorted["ratio"] = df_sorted.apply(compute_ratio, axis=1)
            df_sorted = df_sorted.sort_values(by="ratio", ascending=False, na_position='last')
        
        self.populate_table(df_sorted)
    
    def table_double_click(self, row, column):
        """When a row is double-clicked, show the gene detail dialog."""
        if self.df is None:
            return
        row_data = {}
        col_count = self.table.columnCount()
        headers = [self.table.horizontalHeaderItem(c).text() for c in range(col_count)]
        for c in range(col_count):
            item = self.table.item(row, c)
            row_data[headers[c]] = item.text() if item else ""
        dialog = DetailDialog(row_data, self)
        dialog.exec_()

# --- Main ---

def main():
    preload_file = None
    if len(sys.argv) > 1:
        preload_file = sys.argv[1]
    app = QApplication(sys.argv)
    window = DataExplorer(preload_file)
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()