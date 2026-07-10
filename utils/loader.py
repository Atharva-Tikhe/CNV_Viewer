import pandas as pd
import numpy as np
import os

def load_raw_data(file_path="data/cnv_data.csv"):
    """
    Loads the raw CSV data and renames the first column to 'Source' (values 'DO', 'PO').
    """
    if not os.path.exists(file_path):
        # Fallback to absolute path or other local relative paths
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        file_path = os.path.join(base_dir, "data", "cnv_data.csv")
    
    df = pd.read_csv(file_path)
    
    # The first column is unnamed in the CSV
    if df.columns[0] == "Unnamed: 0" or df.columns[0] == "" or df.columns[0].startswith("Unnamed"):
        df.rename(columns={df.columns[0]: "Source"}, inplace=True)
    else:
        # If it's already named or something else, make sure it is named 'Source'
        df.rename(columns={df.columns[0]: "Source"}, inplace=True)
        
    return df

def get_genes(df):
    """
    Returns the list of genes (columns that are not 'Source' or 'RegId').
    """
    cols = df.columns.tolist()
    if "Source" in cols:
        cols.remove("Source")
    if "RegId" in cols:
        cols.remove("RegId")
    return cols

def get_patients(df):
    """
    Returns the unique RegIds.
    """
    return sorted(df["RegId"].unique().tolist())

def get_long_comparison_df(df):
    """
    Reshapes the data to a long-form comparison table.
    Columns: RegId, Gene, DO_Status, PO_Status, Concordant, Discrepancy_Type
    """
    genes = get_genes(df)
    
    # Split into DO and PO
    do_df = df[df["Source"] == "DO"].copy()
    po_df = df[df["Source"] == "PO"].copy()
    
    # Melt both to long form
    do_long = pd.melt(do_df, id_vars=["RegId"], value_vars=genes, var_name="Gene", value_name="DO_Status")
    po_long = pd.melt(po_df, id_vars=["RegId"], value_vars=genes, var_name="Gene", value_name="PO_Status")

    
    # Merge on RegId and Gene
    merged = pd.merge(do_long, po_long, on=["RegId", "Gene"])
    
    # Determine concordance and discrepancy types
    # DO and PO values are 'N' or 'D'
    # Default to True/False concordance
    merged["Concordant"] = merged["DO_Status"] == merged["PO_Status"]
    
    # Discrepancy types:
    # 'TN' (True Negative): DO='N', PO='N'
    # 'TP' (True Positive): DO='D', PO='D'
    # 'FP' (False Positive): DO='N', PO='D'
    # 'FN' (False Negative): DO='D', PO='N'
    def get_class(row):
        do = row["DO_Status"]
        po = row["PO_Status"]
        if do == "D" and po == "D":
            return "TP"
        elif do == "N" and po == "N":
            return "TN"
        elif do == "N" and po == "D":
            return "FP"
        elif do == "D" and po == "N":
            return "FN"
        return "Unknown"
        
    merged["Class"] = merged.apply(get_class, axis=1)
    
    # Add human readable labels
    class_labels = {
        "TP": "True Positive (Match Deletion)",
        "TN": "True Negative (Match Neutral)",
        "FP": "False Positive (PO Deletion, DO Neutral)",
        "FN": "False Negative (PO Neutral, DO Deletion)"
    }
    merged["Class_Label"] = merged["Class"].map(class_labels)
    
    return merged

def get_patient_wide_comparison(df):
    """
    Returns a wide dataframe where each row is a patient (RegId),
    and columns are the genes showing comparison strings like 'N -> N', 'D -> D', 'N -> D', 'D -> N'.
    """
    long_df = get_long_comparison_df(df)
    
    # Create comparison string
    long_df["Comparison"] = long_df["DO_Status"] + " → " + long_df["PO_Status"]
    
    # Pivot back to wide
    wide_df = long_df.pivot(index="RegId", columns="Gene", values="Comparison").reset_index()
    
    # Add summary columns
    patient_stats = long_df.groupby("RegId").agg(
        Total_Genes=("Gene", "count"),
        Concordant_Genes=("Concordant", "sum"),
        TP=("Class", lambda x: (x == "TP").sum()),
        TN=("Class", lambda x: (x == "TN").sum()),
        FP=("Class", lambda x: (x == "FP").sum()),
        FN=("Class", lambda x: (x == "FN").sum())
    ).reset_index()
    
    patient_stats["Concordance_Rate"] = (patient_stats["Concordant_Genes"] / patient_stats["Total_Genes"] * 100).round(1)
    patient_stats["Has_Discrepancy"] = patient_stats["Concordant_Genes"] < patient_stats["Total_Genes"]
    
    return pd.merge(wide_df, patient_stats, on="RegId")
