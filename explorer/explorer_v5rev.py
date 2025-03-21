import sys
import os
import json
import numpy as np
import pandas as pd
import requests
from PyQt5.QtCore import Qt, QAbstractTableModel, QModelIndex, QThread, pyqtSignal, QPoint
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QTableView, QVBoxLayout, QWidget,
    QFileDialog, QComboBox, QPushButton, QHBoxLayout, QLabel, QDialog, QScrollArea,
    QGroupBox, QMessageBox, QHeaderView, QGridLayout, QTextEdit, QProgressDialog, QLineEdit, QMenu, QDialogButtonBox
)
from PyQt5.QtGui import QColor

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
        self.details_dict = dict(row_data)
        
        main_layout = QVBoxLayout(self)
        
        copy_button = QPushButton("Copy as JSON")
        copy_button.clicked.connect(self.copy_as_json)
        main_layout.addWidget(copy_button)
        
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        main_layout.addWidget(scroll)
        
        content = QWidget()
        scroll.setWidget(content)
        content_layout = QVBoxLayout(content)
        
        exclude_keys = {"fbgn", "arm", "start", "end"}
        basic_fields = [(k, v) for k, v in row_data.items() if k.lower() not in exclude_keys]
        
        if basic_fields:
            grid = QGridLayout()
            grid.setSpacing(10)
            for idx, (key, value) in enumerate(basic_fields):
                widget = self.create_field_box(key, value)
                row = idx // 2
                col = idx % 2
                grid.addWidget(widget, row, col)
            content_layout.addLayout(grid)
        
        fbgn = str(row_data.get("fbgn", "") or "").strip()
        arm = str(row_data.get("arm", "") or "").strip()
        start = row_data.get("start", "")
        end = row_data.get("end", "")
        
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
        
        if ext_layout.count() > 0:
            content_layout.addLayout(ext_layout)
    
    def create_field_box(self, title, content):
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

# --- Location Filter Dialog ---
class LocationFilterDialog(QDialog):
    def __init__(self, arm="", start="", end="", mode="Strict", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Location Filter")
        self.resize(400, 150)
        
        layout = QVBoxLayout(self)
        
        self.armLineEdit = QLineEdit(arm)
        self.armLineEdit.setPlaceholderText("Arm")
        self.startLineEdit = QLineEdit(str(start))
        self.startLineEdit.setPlaceholderText("Start (integer)")
        self.endLineEdit = QLineEdit(str(end))
        self.endLineEdit.setPlaceholderText("End (integer)")
        self.modeCombo = QComboBox()
        self.modeCombo.addItems(["Strict", "Inclusive"])
        if mode in ["Strict", "Inclusive"]:
            self.modeCombo.setCurrentText(mode)
        
        form_layout = QGridLayout()
        form_layout.addWidget(QLabel("Arm:"), 0, 0)
        form_layout.addWidget(self.armLineEdit, 0, 1)
        form_layout.addWidget(QLabel("Start:"), 1, 0)
        form_layout.addWidget(self.startLineEdit, 1, 1)
        form_layout.addWidget(QLabel("End:"), 2, 0)
        form_layout.addWidget(self.endLineEdit, 2, 1)
        form_layout.addWidget(QLabel("Filter Mode:"), 3, 0)
        form_layout.addWidget(self.modeCombo, 3, 1)
        layout.addLayout(form_layout)
        
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Apply | QDialogButtonBox.Cancel)
        self.buttonBox.button(QDialogButtonBox.Apply).clicked.connect(self.apply)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)
    
    def apply(self):
        try:
            self.start_val = int(self.startLineEdit.text().strip())
            self.end_val = int(self.endLineEdit.text().strip())
        except ValueError:
            QMessageBox.warning(self, "Input Error", "Start and End must be integers.")
            return
        self.arm_val = self.armLineEdit.text().strip()
        self.filter_mode = self.modeCombo.currentText()
        self.accept()
    
    def get_filter_params(self):
        return {
            "arm": self.arm_val,
            "start": self.start_val,
            "end": self.end_val,
            "mode": self.filter_mode
        }

# --- Main Data Explorer Window using QTableView ---
class DataExplorer(QMainWindow):
    def __init__(self, preload_file=None):
        super().__init__()
        self.setWindowTitle("Data Explorer")
        self.resize(1000, 600)
        self.df = None
        self.full_df = None  # Store full dataset here
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
        
        # Only Show Location Filter and Clear Filter buttons are included now.
        self.showFilterDialogButton = QPushButton("Show Location Filter")
        self.showFilterDialogButton.clicked.connect(self.open_location_filter_dialog)
        control_layout.addWidget(self.showFilterDialogButton)
        
        self.clearFilterButton = QPushButton("Clear Filter")
        self.clearFilterButton.clicked.connect(self.clear_filter)
        control_layout.addWidget(self.clearFilterButton)
        
        main_layout.addLayout(control_layout)
        self.update_temperature_controls()
        
        self.tableView = QTableView()
        self.tableView.setModel(self.model)
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tableView.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tableView.customContextMenuRequested.connect(self.open_context_menu)
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
        self.full_df = df.copy()  # Save the full dataset
        self.df = df.copy()
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
    
    def open_context_menu(self, pos: QPoint):
        index = self.tableView.indexAt(pos)
        if not index.isValid():
            return
        
        row = index.row()
        menu = QMenu(self)
        detail_action = menu.addAction("Show Detail")
        filter_action = menu.addAction("Filter with Location")
        action = menu.exec_(self.tableView.viewport().mapToGlobal(pos))
        if action == detail_action:
            self.show_detail(row)
        elif action == filter_action:
            self.filter_by_row_location(row)
    
    def show_detail(self, row):
        if self.df is None or row < 0 or row >= len(self.df):
            return
        row_data = self.df.iloc[row].to_dict()
        dialog = DetailDialog(row_data, self)
        dialog.exec_()
    
    def filter_by_row_location(self, row):
        if self.df is None or row < 0 or row >= len(self.df):
            return
        row_data = self.df.iloc[row]
        arm_val = str(row_data.get("arm", "")).strip()
        start_val = row_data.get("start", "")
        end_val = row_data.get("end", "")
        dlg = LocationFilterDialog(arm=arm_val, start=start_val, end=end_val, mode="Strict", parent=self)
        if dlg.exec_():
            params = dlg.get_filter_params()
            self.apply_location_filter(params)
    
    def open_location_filter_dialog(self):
        dlg = LocationFilterDialog(parent=self)
        if dlg.exec_():
            params = dlg.get_filter_params()
            self.apply_location_filter(params)
    
    def apply_location_filter(self, params):
        if self.full_df is None:
            return
        arm_val = params.get("arm", "")
        start_val = params.get("start")
        end_val = params.get("end")
        mode = params.get("mode", "Strict")
        
        filtered_df = self.full_df.copy()
        if arm_val:
            filtered_df = filtered_df[filtered_df['arm'] == arm_val]
        if mode == "Strict":
            filtered_df = filtered_df[(filtered_df['start'] >= start_val) & (filtered_df['end'] <= end_val)]
        else:  # Inclusive mode
            filtered_df = filtered_df[(filtered_df['end'] >= start_val) & (filtered_df['start'] <= end_val)]
        if filtered_df.empty:
            QMessageBox.information(self, "No Results", "No genes found in the specified range.")
            return
        self.df = filtered_df.copy()  # Update current state
        self.model.updateDataFrame(filtered_df)
    
    def clear_filter(self):
        if self.full_df is None:
            return
        self.df = self.full_df.copy()
        self.model.updateDataFrame(self.full_df)
    
def main():
    preload_file = "./output/processed_allele_counts.tsv"
    if len(sys.argv) > 1:
        preload_file = sys.argv[1]
    app = QApplication(sys.argv)
    window = DataExplorer(preload_file)
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()