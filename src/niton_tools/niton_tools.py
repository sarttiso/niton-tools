import numpy as np
import pandas as pd
import geochemdb
import os
from pathlib import Path
import json
from typing import Union

from ipyfilechooser import FileChooser
import ipywidgets as w
from IPython.display import display, clear_output, HTML

from scipy import odr
from scipy import stats

import matplotlib.pyplot as plt

# options
pd.set_option('future.no_silent_downcasting', True)

# load niton_measurement_dict.csv from same directory as this file
current_dir = Path(__file__).parent
meas_dict_df = pd.read_csv(current_dir / 'niton_measurement_dict.csv')
quantities = meas_dict_df['quantity'].unique()

# maps of mean and uncertainty columns for each quantity
quantity_column_map = {}
# maps of measurement units for each quantity
quantity_unit_map = {}
for quantity in quantities:
    cur_df = meas_dict_df[meas_dict_df['quantity'] == quantity]
    quantity_column_map[quantity] = dict(zip(cur_df['type'], cur_df['niton column']))
    quantity_unit_map[quantity] = dict(zip(cur_df['type'], cur_df['unit']))

# helper functions to map element columns to multiindex levels
def build_newcol(row):
    if row["type"] == "mean":
        return (row["quantity"], "standard mean")
    elif row["type"] == "uncertainty":
        return (row["quantity"], "standard 2-Sigma")
    else:
        return (row["quantity"], row["type"])

col_map = {
    row["niton column"]: build_newcol(row)
    for _, row in meas_dict_df.iterrows()
}

def validate_standard_measurements(df: pd.DataFrame):
    """Validate the standard measurements dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to validate.

    Raises
    ------
    ValueError
        If any required columns are missing or if any values are invalid.
    """
    required_columns = ['Reading No', 'Reading Type', 'Time', 'Sample Depth']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
       
    # check that at least one measurement column is present
    measurement_columns = meas_dict_df['niton column'].tolist()
    if not any(col in df.columns for col in measurement_columns):
        raise ValueError("No measurement columns found in dataframe.")
    

def generate_dataframes(standard_file, sheet_name, 
                        geochemdb_path: Union[str, os.PathLike] = None,
                        score_threshold: int = 95):
    """Generate GeochemDB dataframes from a standard run.

    Parameters
    ----------
    standard_file : str
        Path to the standard file.
    sheet_name : str
        Name of the sheet to process.
    geochemdb_path : str, optional
        Path to the GeochemDB SQLite database, by default None. If None, assumes that the database is in the current working directory.
    score_threshold : int, optional
        Score threshold for matching strings in GeochemDB, by default 95.
    """
    df = pd.read_excel(standard_file, sheet_name=sheet_name)
    
    # validate dataframe
    validate_standard_measurements(df)
        
    # connect to geochemdb
    if geochemdb_path is None:
        geochem_path = Path.cwd() / 'standard_database.db'
    else:
        geochem_path = Path(geochemdb_path)
    geodb = geochemdb.GeochemDB(geochem_path)

    # prepare dataframes
    df_aliquots = pd.DataFrame(columns=['aliquot', 'sample', 'material'])
    df_analyses = pd.DataFrame(columns=['analysis', 'aliquot', 'sample', 'date', 'instrument', 'technique'])
    df_measurements = pd.DataFrame(columns=['analysis', 'quantity', 'mean', 'measurement_unit', 'uncertainty', 'uncertainty_unit', 'reference_material'])

    #
    # process aliquots
    #
    # first match aliquots
    aliquot_matches = geodb.matchrows_strings('Aliquots', 
                                              df['Sample Depth'].unique(), 
                                              'aliquot', 
                                              score_threshold=score_threshold)[1]
    # if no matches, error (assume aliquots must already exist in database)
    if len(aliquot_matches) == 0:
        raise ValueError("No matching aliquots found in GeochemDB.")
    # keep only analyses with matched aliquots
    df = df[df['Sample Depth'].isin(aliquot_matches.keys())]
    # populate aliquots dataframe
    df_aliquots['aliquot'] = aliquot_matches.values()
    df_aliquots['sample'] = aliquot_matches.values() # assuming sample name is same as aliquot name
    df_aliquots['material'] = 'powder'  # all standards are currently powders

    #
    # process analyses
    #
    # analysis is yyyy-mm-dd_reading-no
    df_analyses['analysis'] = pd.to_datetime(df['Time']).dt.strftime('%Y-%m-%d') + '_' + df['Reading No'].astype(str)
    # df_analyses['analysis'] = df['Reading No'].astype(str)  
    df_analyses['aliquot'] = df['Sample Depth'].map(aliquot_matches)
    df_analyses['sample'] = df_analyses['aliquot']  # assuming sample name is same as aliquot name
    df_analyses['date'] = pd.to_datetime(df['Time']).astype(str)
    df_analyses['instrument'] = 'Niton XL5 Plus'  # assuming this instrument
    df_analyses['technique'] = df['Reading Type'] # usually ucsb mining

    #
    # process measurements
    #
    # replace "<LOD" with np.nan in niton columns
    df[meas_dict_df['niton column']] = df[meas_dict_df['niton column']].replace(to_replace='<LOD', value=np.nan)

    # iterate over analyses and quantities to populate measurements dataframe
    measurements_list = []
    for _, row in df.iterrows():
        analysis = pd.to_datetime(row['Time']).strftime('%Y-%m-%d') + '_' + str(row['Reading No'])
        for quantity in quantities:
            mean_col = quantity_column_map[quantity].get('mean')
            uncert_col = quantity_column_map[quantity].get('uncertainty')
            # only process if both mean and uncertainty columns exist
            if mean_col and uncert_col:
                mean_value = row[mean_col]
                uncert_value = row[uncert_col]
                # only add measurement if mean_value is not nan
                if pd.notna(mean_value):
                    measurements_list.append({
                        'analysis': analysis,
                        'quantity': quantity,
                        'mean': mean_value,
                        'measurement_unit': quantity_unit_map[quantity]['mean'],
                        'uncertainty': uncert_value,
                        'uncertainty_unit': quantity_unit_map[quantity]['uncertainty'],
                        'reference_material': ''  # no reference material
                    })
    df_measurements = pd.DataFrame(measurements_list)

    return df_measurements, df_analyses, df_aliquots


class StandardUI:
    """Class to handle UI for adding or updating standard measurements to the standard database.
    
    The UI uses ipywidgets and allows the user to locate the standard database as well as the standard measurements spreadsheet.
    Once both files are located and validated, the user can select which sheet to process and then add or update the measurements in the database.
    Prior to adding/updating, the user can verify the matches between aliquots in the spreadsheet and those in the database.

    Attributes
    ----------
    ui : ipywidgets.VBox
        The main UI container widget.
    geodb : geochemdb.GeochemDB
        The GeochemDB object connected to the selected standards database.
    df : pd.DataFrame
        The dataframe containing the selected sheet's data.
    """
    def __init__(self):
        #
        # set up UI elements 
        #

        # database file selection widget (get path to standard database)
        self.standard_db_chooser = FileChooser(start_path='.',  # start in the notebook’s cwd
                                               title='<b>Select standard_database.db</b>',
                                               filter_pattern='standard_database.db',   
                                               show_hidden=False)
        # assign callback function to load and validate the selected database file
        self.standard_db_chooser.register_callback(self.on_standard_db_file_selected)  

        # measurements spreadsheet selection widget (get path to measurements spreadsheet)
        self.standard_meas_chooser = FileChooser(start_path='.',  # start in the notebook’s cwd
                                               title='<b>Select standard_measurements_niton.xlsx</b>',
                                               filter_pattern='standard_measurements_niton.xlsx', 
                                               show_hidden=False)
        
        file_select_box = w.VBox([w.HTML("<b>1. Select Files</b>"),
                                  w.HBox([self.standard_db_chooser, self.standard_meas_chooser],
                                         layout=w.Layout(padding='5px'))],
                                 layout=w.Layout(margin='10px 0px 10px 0px', 
                                                 border='1px solid gray', 
                                                 padding='5px'))
        
        # selection widget for sheet names (populated after spreadsheet is selected)
        self.sheet_select = w.Select(options=[], description='Sheet:', rows=20)
        self.standard_meas_chooser.register_callback(self.on_standard_meas_file_selected)
        # output widget to show dataframe corresponding to selected sheet (only updated after user selects sheet)
        self.sheet_output = w.Output()
        self.sheet_select.observe(self.on_sheet_selected, names='value')
        sheet_select_box = w.VBox([w.HTML("<b>2. Select Sheet to Process</b>"),
                                   w.HBox([self.sheet_select, self.sheet_output], 
                                          layout=w.Layout(padding='5px'))],
                                          layout=w.Layout(margin='10px 0px 10px 0px', 
                                                          border='1px solid gray', 
                                                          padding='5px'))

        # match aliquots with those in the database and show mapping
        self.score_threshold_float = w.BoundedFloatText(value=95, 
                                                        min=0, 
                                                        max=100,
                                                        step=1, 
                                                        description='String Similarity Score:',
                                                        style={'description_width': '130px'},
                                                        layout=w.Layout(width='200px',
                                                                        align_items='flex-start'))
        self.match_button = w.Button(description='Match Standards')
        self.match_output = w.Output()
        self.match_button.on_click(self.on_match_button_clicked)
        aliquot_match_box = w.VBox([w.HTML("<b>3. Match to Standards in Database</b>"),
                                    w.HBox([self.score_threshold_float, 
                                            self.match_button, 
                                            self.match_output],
                                    layout=w.Layout(padding='5px'))],
                                    layout=w.Layout(margin='10px 0px 10px 0px', 
                                                    border='1px solid gray', 
                                                    padding='5px'))

        # button to check if analyses present in database
        self.analysis_check_button = w.Button(description='Check Analyses in DB')
        self.analysis_check_output = w.Output()
        self.analysis_check_button.on_click(self.on_analysis_check_button_clicked)
        analysis_check_box = w.VBox([w.HTML("<b>4. Check Analyses in Database</b>"),
                                     w.HBox([self.analysis_check_button, 
                                             self.analysis_check_output],
                                     layout=w.Layout(padding='5px'))],
                                     layout=w.Layout(margin='10px 0px 10px 0px', 
                                                     border='1px solid gray', 
                                                     padding='5px'))
        
        # section for adding/updating measurements
        self.add_measurements_button = w.Button(description='Add Measurements to DB',
                                                style={'button_color': 'lightgreen'},
                                                layout=w.Layout(width='200px'))
        self.add_measurements_output = w.Output()
        self.add_measurements_button.on_click(self.on_add_measurements_button_clicked)
        add_measurements_box = w.HBox([self.add_measurements_button, 
                                       self.add_measurements_output],
                                      layout=w.Layout(padding='5px'))
        self.update_measurements_button = w.Button(description='Update Measurements in DB',
                                                   style={'button_color': 'lightblue'},
                                                   layout=w.Layout(width='200px'))
        self.update_measurements_output = w.Output()
        self.update_measurements_button.on_click(self.on_update_measurements_button_clicked)
        update_measurements_box = w.HBox([self.update_measurements_button, 
                                          self.update_measurements_output],
                                         layout=w.Layout(padding='5px'))
        add_update_box = w.VBox([w.HTML("<b>5. Add or Update Measurements in Database</b>"),
                                 add_measurements_box,
                                 update_measurements_box],
                                layout=w.Layout(margin='10px 0px 10px 0px', 
                                                border='1px solid gray', 
                                                padding='5px'))

        # combine 
        self.ui = w.VBox([
                    # file selection widgets
                    file_select_box, 
                    # sheet selection and display
                    sheet_select_box,
                    # match aliquots section
                    aliquot_match_box,
                    # check analyses section
                    analysis_check_box,
                    # add/update measurements section
                    add_update_box
                ])
        
        display(self.ui)

    def on_standard_db_file_selected(self, change):
        """Callback for when a standard database file is selected.

        Attempts to create a GeochemDB object with database connection to validate the file.

        Parameters
        ----------
        change : dict
            Dictionary containing information about the change.
        """
        if not self.standard_db_chooser.selected:
            print("No file selected.")
            return
        # path to selected file
        file_path = Path(self.standard_db_chooser.selected)
        if not file_path.exists():
            print(f"File {file_path} does not exist.")
            return
        try:
            # try to connect to the database to validate it
            self.geodb = geochemdb.GeochemDB(file_path)
            print(f"Connected to GeochemDB at {file_path}.")
        except Exception as e:
            print(f"Error connecting to GeochemDB: {e}")

    def on_standard_meas_file_selected(self, change):
        """Callback for when a standard measurements file is selected.

        Parameters
        ----------
        change : dict
            Dictionary containing information about the change.
        """
        if not self.standard_meas_chooser.selected:
            print("No file selected.")
            self.sheet_select.options = []
            self.sheet_select.disabled = True
            return
        # path to selected file
        file_path = Path(self.standard_meas_chooser.selected)
        try:
            # get sheet names
            xls = pd.ExcelFile(file_path)
            self.sheet_select.options = xls.sheet_names
            self.sheet_select.disabled = False
            print(f"Loaded {len(xls.sheet_names)} sheets from {file_path.name}. Please select a sheet to process.")
        except Exception as e:
            print(f"Error loading Excel file: {e}")
            self.sheet_select.options = []
            self.sheet_select.disabled = True

    def on_sheet_selected(self, change):
        """Callback for when a sheet is selected.

        Parameters
        ----------
        change : dict
            Dictionary containing information about the change.
        """
        if not self.standard_meas_chooser.selected or not self.sheet_select.value:
            return
        # path to selected file
        file_path = Path(self.standard_meas_chooser.selected)
        sheet_name = self.sheet_select.value
        try:
            # read the selected sheet into a dataframe
            self.df = pd.read_excel(file_path, sheet_name=sheet_name)
            # validate dataframe
            validate_standard_measurements(self.df)
            # display dataframe in output widget
            with self.sheet_output:
                clear_output()
                display(self.df[['Reading No', 'Reading Type', 'Time', 'Sample Depth']])
            print(f"Loaded sheet '{sheet_name}' with {len(self.df)} rows.")
        except Exception as e:
            with self.sheet_output:
                clear_output()
            print(f"Error loading sheet '{sheet_name}': {e}")

    def on_match_button_clicked(self, b):
        """Callback for when the match button is clicked.

        Parameters
        ----------
        b : Button
            The button that was clicked.
        """
        if not self.standard_db_chooser.selected:
            with self.match_output:
                clear_output()
                print("Please select the standard database file.")
            return
        if not self.standard_meas_chooser.selected or not self.sheet_select.value:
            with self.match_output:
                clear_output()
                print("Please select a standard measurements file and sheet.")
            return
        # path to selected files
        db_path = Path(self.standard_db_chooser.selected)
        meas_path = Path(self.standard_meas_chooser.selected)
        sheet_name = self.sheet_select.value
        score_threshold = self.score_threshold_float.value
        try:
            # load sheet as df
            df = pd.read_excel(meas_path, sheet_name=sheet_name)

            # match aliquots again to show mapping
            aliquot_matches = self.geodb.matchrows_strings('Aliquots', 
                                                      df['Sample Depth'].unique(), 
                                                      'aliquot', 
                                                      score_threshold=score_threshold)[1]
            with self.match_output:
                clear_output()
                if len(aliquot_matches) == 0:
                    print("No matching standards found in GeochemDB.")
                    return
                print("Standard Matches:")
                for orig, matched in aliquot_matches.items():
                    print(f"  {orig} -> {matched}")
                # if aliquots didn't match, print those separately
                unmatched = set(df['Sample Depth'].unique()) - set(aliquot_matches.keys())
                if unmatched:
                    print("\nUnmatched Standards:")
                    for orig in unmatched:
                        print(f"  {orig}")
                # summarize matching; if all matched, say so
                if len(unmatched) == 0:
                    print("\nAll standards matched successfully.")
                
        except Exception as e:
            with self.match_output:
                clear_output()
            print(f"Error processing data: {e}")

    def on_analysis_check_button_clicked(self, b):
        """Callback for when the analysis check button is clicked.

        Parameters
        ----------
        b : Button
            The button that was clicked.
        """
        # if files haven't been selected, inform user
        if not self.standard_db_chooser.selected:
            with self.analysis_check_output:
                clear_output()
                print("Please select the standard database file.")
            return
        if not self.standard_meas_chooser.selected or not self.sheet_select.value:
            with self.analysis_check_output:
                clear_output()
                print("Please select a standard measurements file and sheet.")
            return
        # if sheet hasn't been loaded, inform user
        if self.sheet_select.value is None:
            with self.analysis_check_output:
                clear_output()
                print("Please select a sheet to process.")
            return
        # check that sheet has been loaded into self.df
        if not hasattr(self, 'df'):
            with self.analysis_check_output:
                clear_output()
                print("Please select a sheet to process.")
            return
        # check that geodb has been loaded
        if not hasattr(self, 'geodb'):
            with self.analysis_check_output:
                clear_output()
                print("Please select the standard database file.")
            return
        # search for analyses in database, and print positive and negative matches
        try:
            analyses_str = pd.to_datetime(self.df['Time']).dt.strftime('%Y-%m-%d') + '_' + self.df['Reading No'].astype(str)
            analyses_in_db = self.geodb.matchrows_strings('Analyses', 
                                                          analyses_str, 
                                                          'analysis', 
                                                          score_threshold=100)[1]
            with self.analysis_check_output:
                clear_output()
                if len(analyses_in_db) == 0:
                    print("No analyses from the selected sheet are present in the database. All analyses can be added.")
                    return
                print("Analyses Present in Database:")
                for orig, matched in analyses_in_db.items():
                    print(f"{matched}")
                # if analyses didn't match, print those separately
                unmatched = set(analyses_str) - set(analyses_in_db.keys())
                if unmatched:
                    print("\nAnalyses Not Present in Database (can be added):")
                    for orig in unmatched:
                        print(f"  {orig}")
                # summarize matching; if all matched, say so
                if len(unmatched) == 0:
                    print("\nAll analyses are already present in the database.")
        except Exception as e:
            with self.analysis_check_output:
                clear_output()
            print(f"Error checking analyses: {e}")
    

    def on_add_measurements_button_clicked(self, b):
        """Callback for when the add measurements button is clicked.

        Adds measurements from the selected sheet to the database.

        Parameters
        ----------
        b : Button
            The button that was clicked.
        """
        # if files haven't been selected, inform user
        if not self.standard_db_chooser.selected:
            with self.add_measurements_output:
                clear_output()
                print("Please select the standard database file.")
            return
        if not self.standard_meas_chooser.selected or not self.sheet_select.value:
            with self.add_measurements_output:
                clear_output()
                print("Please select a standard measurements file and sheet.")
            return
        # if sheet hasn't been loaded, inform user
        if self.sheet_select.value is None:
            with self.add_measurements_output:
                clear_output()
                print("Please select a sheet to process.")
            return
        # check that sheet has been loaded into self.df
        if not hasattr(self, 'df'):
            with self.add_measurements_output:
                clear_output()
                print("Please select a sheet to process.")
            return
        # check that geodb has been loaded
        if not hasattr(self, 'geodb'):
            with self.add_measurements_output:
                clear_output()
                print("Please select the standard database file.")
            return
        # path to selected files
        db_path = Path(self.standard_db_chooser.selected)
        meas_path = Path(self.standard_meas_chooser.selected)
        sheet_name = self.sheet_select.value
        score_threshold = self.score_threshold_float.value
        try:
            # generate dataframes from selected sheet
            df_measurements, df_analyses, df_aliquots = generate_dataframes(meas_path, 
                                                                           sheet_name, 
                                                                           geochemdb_path=db_path,
                                                                           score_threshold=score_threshold)
            # add dataframes to database
            with self.add_measurements_output:
                clear_output()
                self.geodb.measurements_add(df_measurements, 
                                            df_analyses, 
                                            df_aliquots, 
                                            score_threshold=score_threshold)

        except Exception as e:
            with self.add_measurements_output:
                clear_output()
            print(f"Error adding measurements: {e}")

    def on_update_measurements_button_clicked(self, b):
        """Callback for when the update measurements button is clicked.

        Updates measurements from the selected sheet in the database.

        Parameters
        ----------
        b : Button
            The button that was clicked.
        """
        # if files haven't been selected, inform user
        if not self.standard_db_chooser.selected:
            with self.update_measurements_output:
                clear_output()
                print("Please select the standard database file.")
            return
        if not self.standard_meas_chooser.selected or not self.sheet_select.value:
            with self.update_measurements_output:
                clear_output()
                print("Please select a standard measurements file and sheet.")
            return
        # if sheet hasn't been loaded, inform user
        if self.sheet_select.value is None:
            with self.update_measurements_output:
                clear_output()
                print("Please select a sheet to process.")
            return
        # check that sheet has been loaded into self.df
        if not hasattr(self, 'df'):
            with self.update_measurements_output:
                clear_output()
                print("Please select a sheet to process.")
            return
        # check that geodb has been loaded
        if not hasattr(self, 'geodb'):
            with self.update_measurements_output:
                clear_output()
                print("Please select the standard database file.")
            return
        # path to selected files
        db_path = Path(self.standard_db_chooser.selected)
        meas_path = Path(self.standard_meas_chooser.selected)
        sheet_name = self.sheet_select.value
        score_threshold = self.score_threshold_float.value
        try:
            # generate dataframes from selected sheet
            df_measurements, _, _ = generate_dataframes(meas_path, 
                                                        sheet_name, 
                                                        geochemdb_path=db_path,
                                                        score_threshold=score_threshold)
            # update measurements in database
            with self.update_measurements_output:
                clear_output()
                self.geodb.measurements_update(df_measurements)

        except Exception as e:
            with self.update_measurements_output:
                clear_output()
            print(f"Error updating measurements: {e}")


class CalibrationEditorUI:
    """UI to create, view, and save calibrations based on standard measurements.
    """
    def __init__(self):
        #
        # set up UI elements
        #
        
        # part 1: file selection
    
        # database file selection widget (get path to standard database)
        self.standard_db_chooser = FileChooser(start_path='.',  # start in the notebook’s cwd
                                               title='<b>Select standard_database.db</b>',
                                               filter_pattern='standard_database.db',   
                                               show_hidden=False)

        # assign callback function to load and validate the selected database file
        self.standard_db_chooser.register_callback(self.on_standard_db_file_selected)

        # select standard definition file
        self.standard_def_chooser = FileChooser(start_path='.',  # start in the notebook’s cwd
                                               title='<b>Select standard_definitions.xlsx</b>',
                                               filter_pattern='standard_definitions.xlsx',
                                               show_hidden=False)
        # assign callback function to load and validate the selected standard definitions file
        self.standard_def_chooser.register_callback(self.on_standard_def_file_selected)

        # assemble file selection widgets
        self.file_chooser_box = w.HBox([self.standard_db_chooser, self.standard_def_chooser],
                                       layout=w.Layout(padding='5px'))
        self.file_chooser_box_container = w.VBox([w.HTML("<b>1. Select Files</b>"),
                                                  self.file_chooser_box],
                                       layout=w.Layout(margin='10px 0px 10px 0px', 
                                                       border='2px solid gray', 
                                                       padding='5px'))

        # part 2: options for filtering standard measurements to include in calibration

        # date range selection                                         
        self.start_date_picker = w.DatePicker(description='Start Date:')
        self.end_date_picker = w.DatePicker(description='End Date:')

        self.date_box = w.VBox([w.HTML('Select Date Range:'),
                                self.start_date_picker, 
                                self.end_date_picker],
                               layout=w.Layout(padding='5px'))
        # create widget for dynamically defined standard checkboxes
        self.standard_box = w.VBox(layout=w.Layout(border='1px solid gray',
                                                  padding='2px',
                                                  width='200px',
                                                  overflow='visible'))  # initially empty
        self.standard_box_container = w.VBox([w.HTML('Select Standards to Include:'),
                                              self.standard_box],
                                              layout=w.Layout(padding='5px'))
        
        # create widget for dynamically defined reading type checkboxes
        self.reading_type_box = w.VBox(layout=w.Layout(border='1px solid gray',
                                                      padding='2px',
                                                      width='200px',
                                                      overflow='visible'))  # initially empty
        self.reading_type_box_container = w.VBox([w.HTML('Select Reading Types to Include:'),
                                                  self.reading_type_box],
                                                  layout=w.Layout(padding='5px'))
        
        # outlier filtering options
        self.outlier_method_dropdown = w.Dropdown(options=['None', 
                                                           'IQR', 
                                                           'Z-Score', 
                                                           'MAD'],
                                                  value='IQR',
                                                  description='Outlier Method:',
                                                  style={'description_width': '200px'})
        self.outlier_threshold_float = w.BoundedFloatText(value=1.5, 
                                                          min=0, 
                                                          max=10,
                                                          step=0.1,
                                                          description='Outlier Threshold:',
                                                          style={'description_width': '200px'},)
        self.outlier_box = w.VBox([w.HTML('Outlier Filtering Options:'),
                                   self.outlier_method_dropdown,
                                   self.outlier_threshold_float],
                                  layout=w.Layout(padding='5px'))
        
        # minimum number of measurements per standard
        self.min_measurements_int = w.BoundedIntText(value=10,
                                                     min=1,
                                                     max=500,
                                                     step=1,
                                                     description='Min Measurements per Standard:',
                                                     style={'description_width': '200px'})
        
        # minimum number of standards per element
        self.min_standards_int = w.BoundedIntText(value=4,
                                                  min=1,
                                                  max=10,
                                                  step=1,
                                                  description='Min Standards per Element:',
                                                  style={'description_width': '200px'})

        # button to run filtering and calibration
        self.run_filter_button = w.Button(description='Filter Standards and Calibrate',
                                         style={'button_color': 'lightgreen'},
                                         layout=w.Layout(width='200px'))
        self.run_filter_button.on_click(self.filter_standards_calibrate)

        
        # assemble filter options
        self.filter_box = w.HBox([w.VBox([self.date_box, 
                                          self.reading_type_box_container,
                                          self.outlier_box,
                                          self.min_measurements_int,
                                          self.min_standards_int],), 
                                  self.standard_box_container,
                                  self.run_filter_button])
        self.filter_box_container = w.VBox([w.HTML("<b>2. Filter Standards for Calibration</b>"),
                                           self.filter_box],
                                           layout=w.Layout(margin='10px 0px 10px 0px', 
                                                           border='2px solid gray', 
                                                           padding='5px'))
        
        
        # part 3: summary standard filtering and calibration information
        
        # create output to provide summary information for each selected standard based on filter values
        self.filter_summary_output = w.Output(layout=w.Layout(border='1px solid gray'))
        self.filter_summary_output_container = w.VBox([w.HTML('Filter Summary:'),
                                                      self.filter_summary_output],
                                                      layout=w.Layout(padding='5px'))
        # output to show calibration summary
        self.calibration_summary_output = w.Output(layout=w.Layout(border='1px solid gray'))
        self.calibration_summary_output_container = w.VBox([w.HTML('Calibration Summary:'),
                                                          self.calibration_summary_output],
                                                          layout=w.Layout(padding='5px'))
        # container for filter summary and calibration summary
        self.summary_container = w.VBox([w.HTML("<b>3. Filter and Calibration Summary</b>"),
                                         w.HBox([self.filter_summary_output_container,
                                                 self.calibration_summary_output_container],)],
                                        layout=w.Layout(margin='10px 0px 10px 0px', 
                                                                 border='2px solid gray', 
                                                                 padding='5px'))

        # part 4: element selection, plotting, calibration

        # widget to hold elements to include in plot
        self.element_box = w.VBox(layout=w.Layout(border='1px solid gray',
                                                    padding='2px',
                                                    width='200px',
                                                    overflow='visible'))  # initially empty
        self.element_box_container = w.VBox([w.HTML('Select Elements to Plot:'),
                                            self.element_box],
                                            layout=w.Layout(padding='5px'))
        
        # set up output to hold plots for selected elements and standard measurements
        self.plotting_output = w.Output()

        # output to show calibration summary
        self.calibration_output = w.Output()

        # container for element selection and plots
        self.plotting_container = w.VBox([w.HTML("<b>4. Plotting</b>"),
                                         w.HBox([self.element_box_container,
                                                 self.plotting_output])],
                                         layout=w.Layout(margin='10px 0px 10px 0px', 
                                                         border='2px solid gray', 
                                                         padding='5px'))
        
        # part 5: calibration saving

        # calibration output folder chooser
        self.calibration_folder_chooser = FileChooser(start_path='.',
                                                      title='<b>Select Calibration Folder</b>',
                                                      show_only_dirs=True)
        
        # calibration file name widget
        self.calibration_name_input = w.Text(description='Calibration File Name:',
                                             style={'button_color': 'lightgreen'},
                                             layout=w.Layout(width='200px'))

        # calibration save button
        self.calibration_save_button = w.Button(description='Save Calibration')
        self.calibration_save_button.on_click(self.calibration_save)

        # output area for feedback
        self.calibration_save_output = w.Output()

        self.calibration_save_container = w.VBox([w.HTML('<b>5. Save Calibration</b>'),
                                                  self.calibration_folder_chooser,
                                                  self.calibration_name_input,
                                                  self.calibration_save_button,
                                                  self.calibration_save_output],
                                                  layout=w.Layout(margin='10px 0px 10px 0px', 
                                                         border='2px solid gray', 
                                                         padding='5px'))

        # assemble UI
        self.ui = w.VBox([
                    # database selection
                    self.file_chooser_box_container,
                    # filter options
                    self.filter_box_container,
                    # filter and calibration summary
                    self.summary_container,
                    # plotting with element selection
                    self.plotting_container,
                    # calibration saving
                    self.calibration_save_container
                ])
        display(self.ui)

    def on_standard_db_file_selected(self, change):
        """Callback for when a standard database file is selected.

        Attempts to create a GeochemDB object with database connection to validate the file.

        Parameters
        ----------
        change : dict
            Dictionary containing information about the change.
        """
        if not self.standard_db_chooser.selected:
            print("No file selected.")
            return
        # path to selected file
        file_path = Path(self.standard_db_chooser.selected)
        if not file_path.exists():
            print(f"File {file_path} does not exist.")
            return
        try:
            # try to connect to the database to validate it
            self.geodb = geochemdb.GeochemDB(file_path)
            # generate standard checkboxes
            self.generate_standard_checkboxes()
            # generate reading type checkboxes
            self.generate_readingtype_checkboxes()
            # set default date range
            self.set_default_date_range()

            print(f"Connected to GeochemDB at {file_path}.")
        except Exception as e:
            print(f"Error connecting to GeochemDB: {e}")

    def on_standard_def_file_selected(self, change):
        """Callback for when a standard definitions file is selected.

        Attempts to load the standard definitions to validate the file.

        Parameters
        ----------
        change : dict
            Dictionary containing information about the change.
        """
        if not self.standard_def_chooser.selected:
            print("No file selected.")
            return
        # path to selected file
        file_path = Path(self.standard_def_chooser.selected)
        if not file_path.exists():
            print(f"File {file_path} does not exist.")
            return
        try:
            # try to load the standard definitions to validate the file
            self.stand_def_df = pd.read_excel(file_path, 
                                              sheet_name='standards',
                                              index_col=0)
            # verify that required columns are present by comparing with niton_measurement_dict
            if not list(self.stand_def_df.columns) >= list(meas_dict_df['niton column']):
                # list missing columns
                missing_cols = set(meas_dict_df['niton column']) - set(self.stand_def_df.columns)
                print(f"Error: Missing required columns in {file_path}: {missing_cols}.")
                return
            # for standards for which elements are not nan, if the corresponding 2-Sigma column is nan, then assign 10% uncertainty
            for ii, el in enumerate(quantities):
                stand_el_idx = ~self.stand_def_df[el].isna()
                cur_sig_col = quantity_column_map[el]['uncertainty']
                stand_el_sig_idx = ~self.stand_def_df[cur_sig_col].isna()
                # indices of standards for current element whose uncertainty must be set to 10%
                idx = stand_el_idx & ~stand_el_sig_idx
                self.stand_def_df.loc[idx, cur_sig_col] = 0.1 * self.stand_def_df.loc[idx, el]

            print(f"Loaded standard definitions from {file_path}.")
        except Exception as e:
            print(f"Error loading standard definitions: {e}")

    def generate_standard_checkboxes(self):
        """Generate checkboxes for selecting standards to include in calibration.
        """
        if not hasattr(self, 'geodb'):
            print("Please select the standard database file.")
            return
        try:
            self.standard_names = self.geodb.get_aliquots()
            self.standard_checkboxes = [w.Checkbox(value=False, 
                                                   indent=False,
                                                   description=name,
                                                   layout=w.Layout(width='100%',
                                                                   height='20px',
                                                                   margin='0px',
                                                                   padding='0px')) for name in self.standard_names]
            
            self.standard_box.children = self.standard_checkboxes
            # set height with buffer, max height 500px
            height_px = min(500, 20 * len(self.standard_checkboxes) + 10) 
            self.standard_box.layout.height = f'{height_px}px'
            # assign callback to filter standard measurements
            for cb in self.standard_checkboxes:
                cb.observe(self.filter_standards, names='value')

        except Exception as e:
            print(f"Error loading standards: {e}")

    def generate_readingtype_checkboxes(self):
        """Generate checkboxes for selecting reading types to include in calibration.
        """
        if not hasattr(self, 'geodb'):
            print("Please select the standard database file.")
            return
        try:
            self.reading_types = pd.read_sql_query('select * from Techniques', self.geodb.con)['name'].tolist()
            self.reading_type_checkboxes = [w.Checkbox(value=False, 
                                                       indent=False,
                                                       description=rt,
                                                       layout=w.Layout(width='100%',
                                                                       height='20px',
                                                                       margin='0px',padding='0px')) for rt in self.reading_types]
            self.reading_type_box.children = self.reading_type_checkboxes
            # set height with buffer, max height 500px
            height_px = min(500, 20 * len(self.reading_type_checkboxes) + 10) 
            self.reading_type_box.layout.height = f'{height_px}px'
            # assign callback to filter standard measurements
            for cb in self.reading_type_checkboxes:
                cb.observe(self.filter_standards, names='value')
        except Exception as e:
            print(f"Error loading reading types: {e}")

    def set_default_date_range(self):
        """Set default date range based on range of dates in self.geodb.

        """
        if not hasattr(self, 'geodb'):
            print("Please select the standard database file.")
            return
        try:
            df_analyses = pd.read_sql_query('select * from Analyses', self.geodb.con)
            df_analyses['date'] = pd.to_datetime(df_analyses['date'], errors='coerce')
            min_date = df_analyses['date'].min()
            max_date = df_analyses['date'].max()
            if pd.isna(min_date) or pd.isna(max_date):
                print("No valid dates found in date column.")
                return
            self.start_date_picker.value = min_date.date()
            self.end_date_picker.value = max_date.date()
            # assign callback to filter standard measurements
            self.start_date_picker.observe(self.filter_standards, names='value')
            self.end_date_picker.observe(self.filter_standards, names='value')
        except Exception as e:
            print(f"Error setting default date range: {e}")

    def set_element_checkboxes(self):
        """Set checkboxes for selecting elements to include in calibration.

        Use elements in self.standard_element_df (i.e., after running filter_standards and process_filtered_standards).
        """
        # check that self.standard_element_df exists
        if not hasattr(self, 'standard_element_df'):
            with self.plotting_output:
                clear_output()
                print("Please filter standards first.")
            return
        
        try:
            # elements present in filtered and processed measurements
            self.elements = self.standard_element_df.index.get_level_values('element').unique().tolist()
            self.elements.sort()
            self.element_checkboxes = [w.Checkbox(value=False,
                                                  indent=False,
                                                  description=el,
                                                  layout=w.Layout(width='100%',
                                                                  height='20px',
                                                                  margin='0px',
                                                                  padding='0px')) for el in self.elements]
            self.element_box.children = self.element_checkboxes

            # assign callback to update plots when elements are selected
            for cb in self.element_checkboxes:
                cb.observe(self.plot_filtered_standards, names='value')
        except Exception as e:
            print(f"Error loading elements: {e}")

    def filter_standards_calibrate(self, change):
        """Callback to filter standards based on selected criteria and perform calibration.

        Criteria include 
        - date range
        - standards to use
        - reading types
        - outlier filtering method and threshold
        - minimum number of measurements per standard
        - minimum number of standards per element

        Parameters
        ----------
        change : dict
            Dictionary containing information about the change.
        """
        if not hasattr(self, 'geodb'):
            with self.filter_summary_output:
                clear_output()
                print("Please select the standard database file.")
            return
        try:
            # get selected date range
            start_date = self.start_date_picker.value
            end_date = self.end_date_picker.value
            if start_date is None or end_date is None:
                with self.filter_summary_output:
                    clear_output()
                    print("Please select a valid start and end date.")
                return
            if start_date > end_date:
                with self.filter_summary_output:
                    clear_output()
                    print("Start date must be before end date.")
                return
            # get list of selected standards
            selected_standards = [cb.description for cb in self.standard_checkboxes if cb.value]
            # ensure at least as many standards are selected as min_standards_int.value
            if len(selected_standards) < self.min_standards_int.value:
                with self.filter_summary_output:
                    clear_output()
                    print(f"Please select at least {self.min_standards_int.value} standards.")
                return
            # get list of selected reading types
            selected_reading_types = [cb.description for cb in self.reading_type_checkboxes if cb.value]
            if len(selected_reading_types) == 0:
                with self.filter_summary_output:
                    clear_output()
                    print("Please select at least one reading type.")
                return
            # query database to get filtered standard measurements
            query = f"""
                SELECT m.*, a.date, a.analysis, al.aliquot, t.name as technique
                FROM Measurements m
                JOIN Analyses a ON m.analysis = a.analysis
                JOIN Aliquots al ON a.aliquot = al.aliquot
                JOIN Techniques t ON a.technique = t.name
                WHERE a.date BETWEEN '{start_date}' AND '{end_date}'
                  AND al.aliquot IN ({','.join(['?']*len(selected_standards))})
                  AND t.name IN ({','.join(['?']*len(selected_reading_types))})
                """
            params = selected_standards + selected_reading_types
            self.filtered_measurements = pd.read_sql_query(query, self.geodb.con, params=params)
            # drop rows with NaN or 0 in "mean" column
            self.filtered_measurements = self.filtered_measurements.dropna(subset=['mean'])
            self.filtered_measurements = self.filtered_measurements[self.filtered_measurements['mean'] != 0]
            # apply outlier filtering if selected
            method = self.outlier_method_dropdown.value
            threshold = self.outlier_threshold_float.value
            self.filtered_measurements = self.filter_outliers(self.filtered_measurements, 
                                                         method=method,
                                                         threshold=threshold)
            
            # group by aliquot and element (multiindex) to summarize number of measurements per element (one columns)
            self.filter_summary = self.filtered_measurements.groupby(['aliquot', 'quantity']).size().unstack(fill_value=0).T
            # create a similar summary dataframe based on the standard definitions to show which elements are defined for each standard
            self.stand_def_summary = self.stand_def_df[self.filter_summary.index].T[self.filter_summary.columns]
            # for each element, create logical array showing which standards have that element defined (not NaN) and with more than self.min_measurements_int.value measurements
            self.standard_filter_df = (self.stand_def_summary.notna()) & (self.filter_summary >= self.min_measurements_int.value)

            # process filtered standards
            self.process_filtered_standards()

            # create element selection checkboxes based on elements present in filtered and processed measurements
            self.set_element_checkboxes()

            # generate calibration
            self.calibrate()

            # display summaries
            self.display_filter_summary()
            self.display_calibration_summary()

            
        except Exception as e:
            print(f"Error filtering and/or calibrating standards: {e}")

    def filter_outliers(self, df, method='IQR', threshold=None):
        """Filter outliers from a dataframe of measurements.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing measurements with columns 'aliquot', 'quantity', and 'mean'.
        method : str, optional
            One of 'IQR', 'Z-Score', 'MAD', or 'None', by default 'IQR'
        threshold : float, optional
            Threshold to use, by default None. If None, then a method-specific threshold is used:
                - 1.5 for IQR
                - 3 for Z-Score
                - 3.5 for MAD
        """
        if method == 'None':
            return df
        filtered_dfs = []
        for (aliquot, quantity), group in df.groupby(['aliquot', 'quantity']):
            if method == 'IQR':
                if threshold is None:
                    threshold = 1.5
                Q1 = group['mean'].quantile(0.25)
                Q3 = group['mean'].quantile(0.75)
                IQR = Q3 - Q1
                lower_bound = Q1 - threshold * IQR
                upper_bound = Q3 + threshold * IQR
                filtered_group = group[(group['mean'] >= lower_bound) & (group['mean'] <= upper_bound)]
            elif method == 'Z-Score':
                if threshold is None:
                    threshold = 3
                mean = group['mean'].mean()
                std = group['mean'].std()
                z_scores = (group['mean'] - mean) / std
                filtered_group = group[abs(z_scores) <= threshold]
            elif method == 'MAD':
                if threshold is None:
                    threshold = 3.5
                median = group['mean'].median()
                mad = np.median(np.abs(group['mean'] - median))
                modified_z_scores = 0.6745 * (group['mean'] - median) / mad
                filtered_group = group[abs(modified_z_scores) <= threshold]
            else:
                raise ValueError(f"Unknown outlier filtering method: {method}")
            filtered_dfs.append(filtered_group)
        if filtered_dfs:
            return pd.concat(filtered_dfs)
        else:
            return pd.DataFrame(columns=df.columns)
        
    def process_filtered_standards(self):
        """Process the filtered standard measurements to compute summary statistics.

        Creates self.standard_element_df with MultiIndex (standard, element) and columns:
            - mean_list: list of mean measurements
            - n: number of measurements
            - sig: standard deviation of measurements
            - max: maximum measurement
            - min: minimum measurement
            - mean: mean of measurements
            - 2.5, 25, 50, 75, 97.5: percentiles of measurements
            - standard mean: known standard value from self.stand_def_df
            - standard 2-Sigma: known standard uncertainty from self.stand_def_df
        """
        if not hasattr(self, 'filtered_measurements'):
            print("Please filter standards first.")
            return
        # aggregate measurements ("mean") by element and aliquot as lists and percentiles, set aliquot and quantities (as 'standard', 'element') as MultiIndex
        agg_funcs = {
            'mean_list': ('mean', list),
            'n': ('mean', 'count'),
            'sig': ('mean', 'std'),
            'max': ('mean', 'max'),
            'min': ('mean', 'min'),
            'mean': ('mean', 'mean'),
        }
        percentiles = [2.5, 25, 50, 75, 97.5]
        for p in percentiles:
            agg_funcs[f'{p}'] = ('mean', lambda x, q=p: np.percentile(x, q))
        self.standard_element_df = self.filtered_measurements.groupby(['aliquot', 'quantity']).agg(**agg_funcs).reset_index()
        self.standard_element_df.rename(columns={'aliquot': 'standard', 'quantity': 'element'}, inplace=True)
        self.standard_element_df.set_index(['standard', 'element'], inplace=True)
        
        # add standard known values and uncertainties from self.stand_def_df
        elem_cols = [c for c in self.stand_def_df.columns if c in col_map]
        stand_def_df_restruct = self.stand_def_df[elem_cols].rename(columns=col_map)
        stand_def_df_restruct.columns = pd.MultiIndex.from_tuples(stand_def_df_restruct.columns)
        stand_def_df_restruct = stand_def_df_restruct.stack(level=0, future_stack=True).rename_axis(['standard', 'element'])
        # perform join
        self.standard_element_df = self.standard_element_df.join(stand_def_df_restruct, how='left')

        # drop any rows with NaN in the standard mean column
        self.standard_element_df = self.standard_element_df.dropna(subset=['standard mean'])
        # filter out standards with fewer than min_measurements_int.value measurements
        self.standard_element_df = self.standard_element_df[self.standard_element_df['n'] >= self.min_measurements_int.value]

        # filter out elements with fewer than min_standards_int.value standards
        self.standard_element_df = self.standard_element_df.groupby('element').filter(lambda x: len(x) >= self.min_standards_int.value)
        self.standard_element_df.index = self.standard_element_df.index.remove_unused_levels()
    
    def display_filter_summary(self):
        """Display summary of filtered standards and elements in the filter_summary_output widget.
        """
        # display summary information for filtered standards
        with self.filter_summary_output:
            clear_output()
            if not hasattr(self, 'filter_summary') or self.filter_summary.empty:
                print("Please filter standards first.")
                return
            if self.filtered_measurements.empty:
                print("No measurements found for the selected criteria.")
                return

            nan_styles = self.stand_def_summary.isna().replace(
                {True: "border: 2px solid red;",
                False: ""}
            )
            # Ensure alignment with summary
            nan_styles = nan_styles.reindex_like(self.filter_summary)
            # also highlight rows with fewer than min_standards_int.value standards in red
            summary_styled = self.filter_summary.style.map(lambda x: 'color: red' if x < self.min_measurements_int.value else 'color: black')\
            .apply(lambda r: [
                # If the row meets the condition → colour the whole row,
                # otherwise leave the cell unchanged (empty string).
                'background-color: mistyrose;' 
                if self.standard_filter_df.loc[r.name].sum() < self.min_standards_int.value
                else ''
                for _ in r   # repeat the same style for every column in this row
                ],
                axis=1        # tell pandas “apply this function to each row”
            )\
            .apply(lambda df: nan_styles, axis=None)
            display(summary_styled)

    def plot_filtered_standards(self, change):
        """Callback for when an element checkbox is changed.

        Updates the plots to show selected elements.

        Parameters
        ----------
        change : dict
            Dictionary containing information about the change.
        """
        # ensure database and standard definitions have been loaded
        if not hasattr(self, 'geodb'):
            with self.plotting_output:
                clear_output()
                print("Please select the standard database file.")
            return
        if not hasattr(self, 'stand_def_df'):
            with self.plotting_output:
                clear_output()
                print("Please select the standard definitions file.")
            return
        # make sure standards have been filtered
        if not hasattr(self, 'filtered_measurements'):
            with self.plotting_output:
                clear_output()
                print("Please filter standards first.")
            return
        # if no measurements after filtering, inform user
        if self.filtered_measurements.empty:
            with self.plotting_output:
                clear_output()
                print("No measurements found for the selected criteria.")
            return
        # ensure calibration has been performed
        if not hasattr(self, 'calibration_df'):
            with self.plotting_output:
                clear_output()
                print("Please perform calibration first.")
            return
        # get list of selected elements
        elements = [cb.description for cb in self.element_checkboxes if cb.value]
        if len(elements) == 0:
            with self.plotting_output:
                clear_output()
                print("Please select at least one element to plot.")
            return
        try:
            num_elements = len(elements)
            fig, ax = plt.subplots(num_elements, 1, figsize=(5, 5*num_elements))
            ax = np.atleast_1d(ax)  # ensure ax is always an array even if num_elements=1
            # create subplots for each selected element
            for ii, el in enumerate(elements):
                # 1:1 line
                ax[ii].axline((0, 0), slope=1, 
                              color='k', linestyle='--',
                              label='1:1')
                # plot calibration with prediction intervals
                cur_x = np.linspace(0, 1.1 * self.standard_element_df.xs(el, level=1)['max'].max(), 100)
                cur_y = self.calibration_df.loc[el]['slope'] * cur_x
                pred_int = prediction_interval(
                    dfdp=[cur_x], # derivative with respect to slope only
                    pcov=self.calibration_df.loc[el]['cov_matrix'],
                    yerr=self.calibration_df.loc[el]['y_err_max'],
                    signif=[68, 95]
                )
                ax[ii].plot(cur_x, cur_y, color='k', label='Calibration fit')
                ax[ii].fill_between(cur_x, cur_y - pred_int[0], cur_y + pred_int[0], color='orange', alpha=0.7, label='68% Prediction Interval')
                ax[ii].fill_between(cur_x, cur_y - pred_int[1], cur_y + pred_int[1], color='moccasin', alpha=0.7, label='95% Prediction Interval')

                # set up violins
                positions = self.stand_def_df.loc[self.standard_element_df.xs(el, level=1).index][el].values
                width = 0.05 * np.max(positions)
                parts = ax[ii].violinplot(self.standard_element_df.xs(el, level=1)['mean_list'].values,
                                positions=positions,
                                widths=width,
                                orientation='horizontal',
                                showmeans=False,
                                showmedians=False,
                                showextrema=False)
                for pc in parts['bodies']:
                    pc.set_alpha(0.6)
                # show box, whiskers
                ax[ii].hlines(positions, 
                              self.standard_element_df.xs(el, level=1)['2.5'].values, 
                              self.standard_element_df.xs(el, level=1)['97.5'].values,
                              color='k', lw=0.5, label='2.5-97.5%')
                ax[ii].hlines(positions, 
                              self.standard_element_df.xs(el, level=1)['25'].values, 
                              self.standard_element_df.xs(el, level=1)['75'].values, 
                              color='k', lw=2, label='25-75%')
                # plot standard uncertainty bars at median value of standard measurements
                medians = self.standard_element_df.xs(el, level=1)['50'].values
                stand_unc = self.stand_def_df.loc[self.standard_element_df.xs(el, level=1).index][f'{el} 2-Sigma'].values
                ax[ii].vlines(medians, positions-stand_unc, positions+stand_unc, color='purple', lw=2, label='Standard uncertainty')
                # label rightmost extents of whiskers with standard names
                for samp, pos, el_max_val in zip(self.standard_element_df.xs(el, level=1).index,
                                                positions, 
                                                self.standard_element_df.xs(el, level=1)['max'].values):
                    ax[ii].text(el_max_val, pos, samp, verticalalignment='center', fontsize=8)

                ax[ii].set_title(f'{el}', fontweight='bold')
                ax[ii].set_ylabel('Known concentration (ppm)')
                ax[ii].set_xlim(0, 1.1 * self.standard_element_df.xs(el, level=1)['max'].max())
                ax[ii].set_ylim(0, 1.1 * np.nanmax(positions+stand_unc))
                ax[ii].grid()

            ax[-1].set_xlabel('Measured concentration (ppm)')
            ax[0].legend()
                
            # display the plots in the output widget
            with self.plotting_output:
                self.plotting_output.clear_output(wait=True)
                display(fig)
                plt.close(fig)  # close the figure to prevent duplicate display in some environments
        except Exception as e:
            print(f"Error updating plots: {e}")

    def calibrate(self):
        """Use the filtered standard measurements to create a calibration.
        
    
        """
        # ensure standards have been selected and filtered
        if not hasattr(self, 'filtered_measurements'):
            with self.calibration_output:
                clear_output()
                print("Please filter standards first.")
            return
        # iterate over all elements and fit line with uncertainty slope and intercept fixed at 0 using measurement uncertainties and standard uncertainties
        calibration_dicts = []
        for el in self.elements:
            # standards used for this element
            cur_el_standards = self.standard_element_df.xs(el, level=1).index.tolist()
            # number of measurements per standard
            cur_el_measurements_per_standard = self.standard_element_df.xs(el, level=1)['n'].values
            # measured
            cur_el_meas_means = self.standard_element_df.xs(el, level=1)['mean'].values
            cur_el_meas_sigs = self.standard_element_df.xs(el, level=1)['sig'].values
            # known
            cur_el_know_means = self.stand_def_df.loc[self.standard_element_df.xs(el, level=1).index][el].values
            cur_el_know_sigs = self.stand_def_df.loc[self.standard_element_df.xs(el, level=1).index][f'{el} 2-Sigma'].values / 2  # convert 2-sigma to 1-sigma
            # skip if no measurements for this element
            if cur_el_meas_means.size == 0:
                continue
            # fit line with slope only, intercept fixed at 0 using orthogonal distance regression
            try:
                # create a Model for ODR
                linear_model = odr.Model(linear_func)

                # create a RealData object using the measured and known values along with their standard deviations
                data = odr.RealData(cur_el_meas_means, 
                                    cur_el_know_means, 
                                    sx=cur_el_meas_sigs, 
                                    sy=cur_el_know_sigs)

                # set up ODR with the model and data
                odr_instance = odr.ODR(data, linear_model, beta0=[1.0])  # initial guess for slope

                # run the regression
                odr_output = odr_instance.run()

                # extract slope and its standard error
                slope = odr_output.beta[0]
                slope_unc = odr_output.sd_beta[0]
                y_err_max = np.max(np.abs(odr_output.eps))

                # reduced variance 95% limit (chi2 level for given dof at 95% confidence level)
                dof = len(cur_el_meas_means)-1
                red_var_95 = stats.chi2.ppf(0.95, dof)/dof

                # record results
                calibration_dicts.append({'element': el,
                                           'slope': slope,
                                           'slope_unc': slope_unc,
                                           'y_err_max': y_err_max,
                                           'reduced variance': odr_output.res_var,
                                           'reduced variance 95%': red_var_95,
                                           'cov_matrix': odr_output.cov_beta,
                                           'standards_used': cur_el_standards,
                                           'measurements_per_standard': cur_el_measurements_per_standard})
            except Exception as e:
                print(f"Error fitting calibration for element {el}: {e}")
                continue

        # create dataframe from calibration_dicts
        self.calibration_df = pd.DataFrame(calibration_dicts).set_index('element')

    def display_calibration_summary(self):
        """Display summary of calibration in the calibration_summary_output widget.
        """
        with self.calibration_summary_output:
            clear_output()
            if not hasattr(self, 'calibration_df') or self.calibration_df.empty:
                print("No calibration available. Please filter standards and generate calibration first.")
                return
            calibration_summary = self.calibration_df[['slope', 
                                                          'slope_unc', 
                                                          'reduced variance',
                                                          'reduced variance 95%']]
                
            # Helper function to style a row based on a condition between two columns
            def style_variance_row(row):
                # Default is no style for any cell in this row
                styles = [''] * len(row)
                # Apply style to 'reduced variance' cell if condition is met
                if row['reduced variance'] > row['reduced variance 95%']:
                    col_idx = row.index.get_loc('reduced variance')
                    styles[col_idx] = 'background-color: mistyrose'
                return styles

            calibration_summary_styled = calibration_summary.style.map(
                # 1. This first map for single-column styling was correct and remains.
                #    It styles 'slope' cells based only on their own value.
                lambda x: 'background-color: mistyrose' if abs(x - 1) > 0.3 else '', 
                subset=['slope']
            ).apply(
                # 2. This apply is the new part for multi-column styling.
                #    It styles 'reduced variance' cells based on another column.
                style_variance_row, 
                axis=1
            )
                                                        
            display(calibration_summary_styled)
    
    def calibration_save(self, change):
        """Save calibration as json file to use for correction unknown standards.
        """
        # ensure calibration exists
        if not hasattr(self, 'calibration_df') or self.calibration_df.empty:
            with self.calibration_save_output:
                clear_output()
                print("No calibration available. Please filter standards and generate calibration first.")
            return
        # create dictionary with standard filtering values
        filter_metadata_dict = {
            'start_date': self.start_date_picker.value.isoformat(),
            'end_date': self.end_date_picker.value.isoformat(),
            'outlier_method': self.outlier_method_dropdown.value,
            'outlier_threshold': self.outlier_threshold_float.value,
            'min_measurements_per_standard': self.min_measurements_int.value,
            'min_standards_per_element': self.min_standards_int.value,
            'reading_types': [cb.description for cb in self.reading_type_checkboxes if cb.value],
            'standards': [cb.description for cb in self.standard_checkboxes if cb.value],
        }
        # dictionary form of filter_summary
        filter_summary_dict = self.filter_summary.to_dict()
        # dictionary form of calibration_df
        calibration_df_lists = self.calibration_df.map(lambda x: x.tolist() if isinstance(x, np.ndarray) else x)
        calibration_dict = calibration_df_lists.to_dict()

        # final dict
        calibration_output_dict = {
            'standard_filter_metadata': filter_metadata_dict,
            'standard_filter_summary': filter_summary_dict,
            'calibration': calibration_dict
        }

        # save as json
        file_name = self.calibration_name_input.value
        if not file_name.endswith('.json'):
            file_name += '.json'
        file_path = Path(self.calibration_folder_chooser.selected) / file_name

        try:
            with open(file_path, 'w') as f:
                json.dump(calibration_output_dict, f, indent=4)
                with self.calibration_save_output:
                    clear_output()
                    print('Calibration saved!')
        except Exception as e:
            with self.calibration_save_output:
                clear_output()
                print(f"Error saving calibration: {e}")


def prediction_interval(dfdp, pcov, yerr, signif=95):
    """Prediction interval for ODR fits.

    Taken from https://stackoverflow.com/questions/60889680/how-to-plot-1-sigma-prediction-interval-for-scipy-odr

    Parameters
    ----------
    dfdp : array-like
        Derivatives of the fit function with respect to the fit parameters at the new x value(s).
        Each element of the list/array corresponds to one parameter and contains an array of derivatives at each x value.
    pcov : 2D array-like
        Covariance matrix of the fit parameters (ODR_output.cov_beta).
    yerr : float or array-like
        Maximum y error of the data points (np.max(ODR_output.eps)).
    signif : float or array-like, optional
        Significance level (e.g., 95 for 95% prediction interval), by default 95.
        If array-like, then output is calculated for each significance level.
    
    Returns
    -------
    float or array-like
        Prediction band offset for a new measurement at the given x value(s). Shape is (n,) if signif is a float, or (len(signif), n) if signif is array-like.
    """
    n_params = len(dfdp)  # number of fit parameters
    n = len(dfdp[0])  # number of data points

    # convert provided value to a significance level (e.g. 95. -> 0.05), then calculate alpha
    signif = np.atleast_1d(signif)
    alpha = 1. - (1 - signif / 100.0) / 2

    # student’s t test
    tval = stats.t.ppf(alpha, n - n_params)

    # process covariance matrix and derivatives
    d = np.zeros(n)
    for j in range(n_params):
        for k in range(n_params):
            d += dfdp[j] * dfdp[k] * pcov[j,k]

    # return prediction band offset for a new measurement
    return np.squeeze(tval.reshape(-1, 1) * np.sqrt(yerr**2 + d))

# define linear function with intercept fixed at 0
def linear_func(B, x):
    return B[0] * x  # B[0] is the slope