import re
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

def change_sample_id_format(sample_id, mb=True):
    """
    params: sample_id:  str
                        Sample ID. eg:- A01-33A (Or A01.33A in case of Poultry samples)
            mb :        boolean {True, False}
                        Whether the data is microbiome or not
    returns: new_index : str
                        Format A1-33A
    """
    sample_id = str(sample_id)
    if mb:
        try:
            farmflock, sample_alpha_num = sample_id.split(".") # A01.33A -> A1 and 33A
        except:
            return
    else:
        farmflock, sample_alpha_num = sample_id.split("-") # A01-33A -> A1 and 33A
    
    farm, flock, _ = re.split('(\d+)',farmflock) # A01 -> A and 01
    flock = str(int(flock)) # 01 to 1

    try:
        sample_num = str(int(sample_alpha_num)) 
        new_index = farm+flock+"-"+sample_num
    except:
        sample_num = str(int(sample_alpha_num[:-1]))
        sample_alpha = sample_alpha_num[-1]
        new_index = farm+flock+"-"+sample_num+sample_alpha
    return new_index

def sheet_transpose(df_ref):
    df_new = df_ref.transpose()
    df_new.index.name = "SampleID"
    return df_new

def get_bacterias_and_samples_list(sheets_list):
    # Returns the list of unique bacterias and unified sample list
    all_indices = []
    all_bacterias = []
    for sheet in sheets_list:
        sheet_index = sheet.index.tolist()
        sheet_cols = sheet.columns.tolist()
        for sample in sheet_index:
            if sample not in all_indices:
                all_indices.append(sample)
        for bacteria in sheet_cols:
            if bacteria not in all_bacterias:
                all_bacterias.append(bacteria)
    
    return all_bacterias, all_indices

def unify_mb_data(filename, num_sheets=8):
    # Get unique bacteria and index names
    sheets_list = []
    for i in range(num_sheets):
        df = pd.read_excel(filename, sheet_name=i, skiprows=1, index_col="#OTU ID")
        df = sheet_transpose(df)
        sheets_list.append(df)

    col_names, indexes = get_bacterias_and_samples_list(sheets_list)

    # New dataframe
    new_df = pd.DataFrame(columns=col_names, index=indexes)
    new_df.index.name = "SampleID"

    # Reassigning
    for i in range(num_sheets):
        print("Sheet number: ", i)
        df = pd.read_excel(filename, sheet_name=i, skiprows=1, index_col="#OTU ID")
        for col in df.columns:
            for row in df.index:
                new_df[row][col] = df[col][row]

    new_df = new_df.fillna(0)
    # Index re formatting
    indices_new_df = new_df.index.tolist()
    indices_formatted = [change_sample_id_format(x) for x in indices_new_df]
    new_df.index = indices_formatted
    return new_df
    
def standard_scale(df_ref):
    """
    Standar scales and zscales data with the column axis
    params: df_ref: pandas dataframe
                    Microbiome data that is unified
    """
    
    def scale2range(col):
        def scale(x):
            scaled = (x-np.min(col))/(np.max(col)-np.min(col)) if x is not None else np.nan
            return scaled

        scaled = [scale(x) for x in col]
        return scaled

    scaler = StandardScaler()
    cols = df_ref.columns
    samples = df_ref.index
    
    # New dataframe
    df_ret = pd.DataFrame(columns=cols, index=samples)
    
    for col in cols:
        data = np.array(df_ref[col].tolist())
        data = data.reshape((-1, 1))
        scaler.fit(data)
        converted = scaler.transform(data)
        converted = np.squeeze(converted)
        converted = scale2range(converted)
        df_ret[col] = converted
    return df_ret
