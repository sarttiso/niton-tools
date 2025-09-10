import numpy as np
import pandas as pd
import geochemdb
import os
from pathlib import Path
from typing import Union

from ipyfilechooser import FileChooser
import ipywidgets as w
from IPython.display import display, clear_output

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
    df_analyses['analysis'] = df['Reading No'].astype(str)
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
        analysis = str(row['Reading No'])
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
                    print("No matching aliquots found in GeochemDB.")
                    return
                print("Aliquot Matches:")
                for orig, matched in aliquot_matches.items():
                    print(f"  {orig} -> {matched}")
                # if aliquots didn't match, print those separately
                unmatched = set(df['Sample Depth'].unique()) - set(aliquot_matches.keys())
                if unmatched:
                    print("\nUnmatched Aliquots:")
                    for orig in unmatched:
                        print(f"  {orig}")
                # summarize matching; if all matched, say so
                if len(unmatched) == 0:
                    print("\nAll aliquots matched successfully.")
                
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
            analyses_in_db = self.geodb.matchrows_strings('Analyses', 
                                                          self.df['Reading No'].astype(str).unique(), 
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
                unmatched = set(self.df['Reading No'].astype(str).unique()) - set(analyses_in_db.keys())
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