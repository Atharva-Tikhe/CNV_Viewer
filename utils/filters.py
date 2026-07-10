import pandas as pd

def apply_filters(long_df, selected_genes=None, selected_patients=None, concordance_filter="All", gene_deletion_filter=None):
    """
    Applies filters to the long-form dataframe.
    
    concordance_filter: "All", "Concordant Only", "Discrepant Only"
    gene_deletion_filter: dict of {gene_name: [list of 'DO_D', 'DO_N', 'PO_D', 'PO_N']} to enforce specific patterns
    """
    df = long_df.copy()
    
    # 1. Filter by Genes
    if selected_genes:
        df = df[df["Gene"].isin(selected_genes)]
        
    # 2. Filter by Patients (RegId)
    if selected_patients:
        df = df[df["RegId"].isin(selected_patients)]
        
    # 3. Filter by Concordance
    if concordance_filter == "Concordant Only":
        df = df[df["Concordant"] == True]
    elif concordance_filter == "Discrepant Only":
        df = df[df["Concordant"] == False]
        
    # 4. Filter by specific gene deletion conditions
    # e.g., if we want to filter patients who have specific deletion patterns.
    # Since this requires patient-level filtering (which affects all genes for that patient),
    # we first find the matching RegIds and then filter the dataframe.
    if gene_deletion_filter:
        matching_patients = set(df["RegId"].unique())
        
        for gene, conditions in gene_deletion_filter.items():
            if not conditions:
                continue
            
            # Find patient IDs meeting this gene's conditions
            gene_df = long_df[long_df["Gene"] == gene]
            gene_matches = set()
            
            for cond in conditions:
                if cond == "DO Deletion (D)":
                    cond_matches = gene_df[gene_df["DO_Status"] == "D"]["RegId"].unique()
                elif cond == "DO Neutral (N)":
                    cond_matches = gene_df[gene_df["DO_Status"] == "N"]["RegId"].unique()
                elif cond == "PO Deletion (D)":
                    cond_matches = gene_df[gene_df["PO_Status"] == "D"]["RegId"].unique()
                elif cond == "PO Neutral (N)":
                    cond_matches = gene_df[gene_df["PO_Status"] == "N"]["RegId"].unique()
                else:
                    cond_matches = []
                
                gene_matches.update(cond_matches)
            
            matching_patients.intersection_update(gene_matches)
            
        df = df[df["RegId"].isin(matching_patients)]
        
    return df
