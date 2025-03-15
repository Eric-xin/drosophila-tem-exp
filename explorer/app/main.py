import sys
import os
import json
import numpy as np
import pandas as pd
import requests
from PyQt5.QtCore import Qt, QAbstractTableModel, QModelIndex, QThread, pyqtSignal
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QTableView, QVBoxLayout, QWidget,
    QFileDialog, QComboBox, QPushButton, QHBoxLayout, QLabel, QDialog, QScrollArea,
    QGroupBox, QMessageBox, QHeaderView, QGridLayout, QTextEdit, QProgressDialog
)
from PyQt5.QtGui import QColor

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for Py2app """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS.
        base_path = sys._MEIPASS  # for PyInstaller, similar for py2app you'll use sys.argv[0]
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

# Use resource_path to set the preload file:
DEFAULT_TSV = resource_path("data/data.tsv")

# --- Worker Thread Class ---
class Worker(QThread):
    result_signal = pyqtSignal(object)
    error_signal = pyqtSignal(str)
    
    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
    
    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
            self.result_signal.emit(result)
        except Exception as e:
            self.error_signal.emit(str(e))

# --- Pandas Model for QTableView ---
class PandasModel(QAbstractTableModel):
    def __init__(self, df=pd.DataFrame(), parent=None):
        super().__init__(parent)
        self._df = df
        
    def rowCount(self, parent=QModelIndex()):
        return len(self._df)
    
    def columnCount(self, parent=QModelIndex()):
        return len(self._df.columns)
    
    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        value = self._df.iloc[index.row(), index.column()]
        if role == Qt.DisplayRole:
            if pd.isna(value):
                return "N/A"
            else:
                col = self._df.columns[index.column()]
                if col in ['start', 'end']:
                    return str(int(value))
                if isinstance(value, float):
                    return f"{value:.4f}"
                return str(value)
        if role == Qt.TextAlignmentRole:
            if isinstance(value, (int, float)) and not pd.isna(value):
                return Qt.AlignRight | Qt.AlignVCenter
            else:
                return Qt.AlignLeft | Qt.AlignVCenter
        if role == Qt.ForegroundRole and pd.isna(value):
            return QColor('gray')
        return None
    
    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            try:
                return self._df.columns[section]
            except IndexError:
                return None
        else:
            return str(section)
    
    def updateDataFrame(self, df):
        self.beginResetModel()
        self._df = df
        self.endResetModel()

# --- FlyBase API functions with improved error handling ---
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

# --- Utility function for formatting sequence with colors ---
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

# --- Detail Dialog using QTextEdit for multi-line display ---
class DetailDialog(QDialog):
    def __init__(self, row_data, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Gene Detail")
        self.resize(700, 500)
        # Copy the passed-in data
        self.details_dict = dict(row_data)
        
        # Main vertical layout for the dialog
        main_layout = QVBoxLayout(self)
        
        # "Copy as JSON" button at the top
        copy_button = QPushButton("Copy as JSON")
        copy_button.clicked.connect(self.copy_as_json)
        main_layout.addWidget(copy_button)
        
        # Scroll area to hold all content
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        main_layout.addWidget(scroll)
        
        # Container widget for the scroll area
        content = QWidget()
        scroll.setWidget(content)
        content_layout = QVBoxLayout(content)
        
        # -------------------------
        # First: Basic row data in grid (2 boxes per row)
        # Exclude keys that we plan to show as external details.
        exclude_keys = {"fbgn", "arm", "start", "end"}
        basic_fields = [(k, v) for k, v in row_data.items() if k.lower() not in exclude_keys]
        
        if basic_fields:
            grid = QGridLayout()
            grid.setSpacing(10)
            # Add each basic field to the grid (2 columns per row)
            for idx, (key, value) in enumerate(basic_fields):
                widget = self.create_field_box(key, value)
                row = idx // 2
                col = idx % 2
                grid.addWidget(widget, row, col)
            content_layout.addLayout(grid)
        
        # -------------------------
        # Second: External details (Gene Summary and Gene Sequence) in vertical layout (one per row)
        # Convert fbgn and location-related values to string safely.
        fbgn = str(row_data.get("fbgn", "") or "").strip()
        arm = str(row_data.get("arm", "") or "").strip()
        start = row_data.get("start", "")
        end = row_data.get("end", "")
        
        # Create a vertical layout for external details
        ext_layout = QVBoxLayout()
        if fbgn and fbgn.startswith("FBgn"):
            summary = fetch_id_abt(fbgn)
            self.details_dict["Gene Summary"] = summary
            summary_box = self.create_field_box("Gene Summary", summary)
            ext_layout.addWidget(summary_box)
            
            seq = fetch_sequence_id(fbgn)
            self.details_dict["Gene Sequence"] = seq
            seq_box = self.create_sequence_box("Gene Sequence", seq)
            ext_layout.addWidget(seq_box)
        elif arm and start and end:
            location = f"{arm}:{start}..{end}"
            seq = fetch_sequence_location("dmel", location)
            self.details_dict["Gene Sequence"] = seq
            seq_box = self.create_sequence_box(f"Gene Sequence (from {location})", seq)
            ext_layout.addWidget(seq_box)
        
        # If external details were added, put a separator (optional) and add them to main layout
        if ext_layout.count() > 0:
            content_layout.addLayout(ext_layout)
    
    def create_field_box(self, title, content):
        """Create a styled box (using QTextEdit) for multi-line display of a field."""
        box = QGroupBox(title)
        layout = QVBoxLayout(box)
        text_edit = QTextEdit()
        text_edit.setPlainText(str(content))
        text_edit.setReadOnly(True)
        text_edit.setFrameStyle(0)
        text_edit.setLineWrapMode(QTextEdit.WidgetWidth)
        layout.addWidget(text_edit)
        return box
    
    def create_sequence_box(self, title, sequence):
        """Create a box for gene sequence with ATCG coloring, using QTextEdit."""
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
        json_str = json.dumps(self.details_dict, indent=4)
        clipboard = QApplication.clipboard()
        clipboard.setText(json_str)
        QMessageBox.information(self, "Copied", "Data copied as JSON to clipboard.")

# --- Main Data Explorer Window using QTableView ---
class DataExplorer(QMainWindow):
    def __init__(self, preload_file=None):
        super().__init__()
        self.setWindowTitle("Data Explorer")
        self.resize(1000, 600)
        self.df = None
        self.model = PandasModel(pd.DataFrame())
        
        main_widget = QWidget(self)
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout(main_widget)
        
        control_layout = QHBoxLayout()
        self.loadButton = QPushButton("Load TSV Data")
        self.loadButton.clicked.connect(self.open_file_dialog)
        control_layout.addWidget(self.loadButton)
        
        self.sortLabel = QLabel("Sort Mode:")
        control_layout.addWidget(self.sortLabel)
        
        self.sortModeCombo = QComboBox()
        self.sortModeCombo.addItems([
            "Sort by Location",
            "Sort by Temperature Reading",
            "Sort by Temperature Ratio"
        ])
        self.sortModeCombo.currentIndexChanged.connect(self.on_sort_mode_change)
        control_layout.addWidget(self.sortModeCombo)
        
        self.tempSortLabel = QLabel("Temperature Reading:")
        self.tempSortCombo = QComboBox()
        self.tempSortCombo.addItems(["temp13_avg", "temp18_avg", "temp23_avg", "temp29_avg"])
        control_layout.addWidget(self.tempSortLabel)
        control_layout.addWidget(self.tempSortCombo)
        
        self.ratioNumLabel = QLabel("Temperature Numerator:")
        self.ratioNumCombo = QComboBox()
        self.ratioNumCombo.addItems(["temp13_avg", "temp18_avg", "temp23_avg", "temp29_avg"])
        self.ratioDenomLabel = QLabel("Temperature Denom:")
        self.ratioDenomCombo = QComboBox()
        self.ratioDenomCombo.addItems(["temp13_avg", "temp18_avg", "temp23_avg", "temp29_avg"])
        control_layout.addWidget(self.ratioNumLabel)
        control_layout.addWidget(self.ratioNumCombo)
        control_layout.addWidget(self.ratioDenomLabel)
        control_layout.addWidget(self.ratioDenomCombo)
        
        self.sortButton = QPushButton("Sort Data")
        self.sortButton.clicked.connect(self.sort_data)
        control_layout.addWidget(self.sortButton)
        
        main_layout.addLayout(control_layout)
        self.update_temperature_controls()
        
        self.tableView = QTableView()
        self.tableView.setModel(self.model)
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tableView.doubleClicked.connect(self.table_double_click)
        main_layout.addWidget(self.tableView)
        
        if preload_file is not None:
            self.load_file(preload_file)
    
    def update_temperature_controls(self):
        mode = self.sortModeCombo.currentText()
        if mode == "Sort by Temperature Ratio":
            self.tempSortLabel.setVisible(False)
            self.tempSortCombo.setVisible(False)
            self.ratioNumLabel.setVisible(True)
            self.ratioNumCombo.setVisible(True)
            self.ratioDenomLabel.setVisible(True)
            self.ratioDenomCombo.setVisible(True)
        elif mode == "Sort by Temperature Reading":
            self.tempSortLabel.setVisible(True)
            self.tempSortCombo.setVisible(True)
            self.ratioNumLabel.setVisible(False)
            self.ratioNumCombo.setVisible(False)
            self.ratioDenomLabel.setVisible(False)
            self.ratioDenomCombo.setVisible(False)
        else:
            self.tempSortLabel.setVisible(False)
            self.tempSortCombo.setVisible(False)
            self.ratioNumLabel.setVisible(False)
            self.ratioNumCombo.setVisible(False)
            self.ratioDenomLabel.setVisible(False)
            self.ratioDenomCombo.setVisible(False)
    
    def on_sort_mode_change(self, index):
        self.update_temperature_controls()
    
    def open_file_dialog(self):
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open TSV File", "", "TSV Files (*.tsv *.txt);;All Files (*)", options=options
        )
        if filename:
            self.load_file(filename)
    
    def load_file(self, filename):
        def load():
            df = pd.read_csv(filename, sep="\t")
            for col in ['start', 'end']:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
            for col in ['temp13_avg', 'temp18_avg', 'temp23_avg', 'temp29_avg']:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            return df
        
        self.worker = Worker(load)
        self.worker.result_signal.connect(self.on_file_loaded)
        self.worker.error_signal.connect(lambda err: QMessageBox.critical(self, "Error", f"Error loading file: {err}"))
        self.worker.start()
    
    def on_file_loaded(self, df):
        self.df = df
        self.model.updateDataFrame(df)
    
    def sort_data(self):
        if self.df is None:
            return
        mode = self.sortModeCombo.currentText()
        df_copy = self.df.copy()
        def do_sort():
            if mode == "Sort by Location":
                return df_copy.sort_values(by=["arm", "start", "end"], ascending=True)
            elif mode == "Sort by Temperature Reading":
                sort_col = self.tempSortCombo.currentText()
                return df_copy.sort_values(by=sort_col, ascending=False, na_position='last')
            elif mode == "Sort by Temperature Ratio":
                num_col = self.ratioNumCombo.currentText()
                denom_col = self.ratioDenomCombo.currentText()
                if num_col == denom_col:
                    raise Exception("Please select two different temperature columns for the ratio.")
                with pd.option_context('mode.use_inf_as_na', True):
                    df_copy["ratio"] = df_copy[num_col] / df_copy[denom_col].replace(0, pd.NA)
                return df_copy.sort_values(by="ratio", ascending=False, na_position='last')
            else:
                return df_copy
        
        self.progress_dialog = QProgressDialog("Sorting data...", None, 0, 0, self)
        self.progress_dialog.setWindowTitle("Sorting")
        self.progress_dialog.setCancelButton(None)
        self.progress_dialog.setModal(True)
        self.progress_dialog.show()
        
        self.sort_worker = Worker(do_sort)
        self.sort_worker.result_signal.connect(self.on_sort_finished)
        self.sort_worker.error_signal.connect(lambda err: QMessageBox.critical(self, "Error", f"Error sorting data: {err}"))
        self.sort_worker.start()
    
    def on_sort_finished(self, df_sorted):
        self.progress_dialog.close()
        self.df = df_sorted
        self.model.updateDataFrame(df_sorted)
    
    def table_double_click(self, index):
        if not index.isValid():
            return
        row = index.row()
        row_data = self.df.iloc[row].to_dict()
        dialog = DetailDialog(row_data, self)
        dialog.exec_()

def main():
    preload_file = DEFAULT_TSV
    if len(sys.argv) > 1:
        preload_file = sys.argv[1]
    app = QApplication(sys.argv)
    window = DataExplorer(preload_file)
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()