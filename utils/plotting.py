import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np

# Set standard colors for consistency across the application
COLOR_MAP = {
    "TP": "#10ac84",      # Teal/Green - Match Deletion
    "TN": "#dfe6e9",      # Light Grey - Match Neutral
    "FP": "#ee5253",      # Soft Red - False Positive
    "FN": "#ff9f43"       # Warm Orange - False Negative
}

COLOR_LABELS = {
    "TP": "True Positive (Match Deletion)",
    "TN": "True Negative (Match Neutral)",
    "FP": "False Positive (PO Deletion, DO Neutral)",
    "FN": "False Negative (PO Neutral, DO Deletion)"
}

def plot_deletion_frequency(long_df):
    """
    Creates a grouped bar chart showing deletion frequency ('D' status) per gene
    for Diagnostic Output (DO) and Pipeline Output (PO).
    """
    # Filter for deletions
    do_dels = long_df[long_df["DO_Status"] == "D"].groupby("Gene").size().reset_index(name="DO")
    po_dels = long_df[long_df["PO_Status"] == "D"].groupby("Gene").size().reset_index(name="PO")
    
    # Get all genes to ensure we have a complete list
    all_genes = sorted(long_df["Gene"].unique())
    freq_df = pd.DataFrame({"Gene": all_genes})
    
    # Merge
    freq_df = freq_df.merge(do_dels, on="Gene", how="left").fillna(0)
    freq_df = freq_df.merge(po_dels, on="Gene", how="left").fillna(0)
    
    # Calculate percentage
    num_patients = long_df["RegId"].nunique()
    freq_df["DO_Pct"] = (freq_df["DO"] / num_patients * 100).round(1)
    freq_df["PO_Pct"] = (freq_df["PO"] / num_patients * 100).round(1)
    
    # Melt for plotly
    plot_df = pd.melt(
        freq_df, 
        id_vars=["Gene"], 
        value_vars=["DO_Pct", "PO_Pct"], 
        var_name="Output Type", 
        value_name="Deletion Frequency (%)"
    )
    
    plot_df["Output Type"] = plot_df["Output Type"].map({
        "DO_Pct": "Diagnostic (DO)",
        "PO_Pct": "Pipeline (PO)"
    })
    
    fig = px.bar(
        plot_df,
        x="Gene",
        y="Deletion Frequency (%)",
        color="Output Type",
        barmode="group",
        color_discrete_map={
            "Diagnostic (DO)": "#34495e",  # Dark slate
            "Pipeline (PO)": "#3498db"     # Premium blue
        },
        title="Deletion Frequency (%) per Gene: DO vs. PO",
        text="Deletion Frequency (%)"
    )
    
    fig.update_traces(textposition='outside', textfont_size=10)
    fig.update_layout(
        hovermode="x unified",
        yaxis_title="Deletion Frequency (%)",
        xaxis_title="Gene",
        legend_title="Source",
        margin=dict(l=40, r=40, t=60, b=40),
        height=450,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        yaxis=dict(gridcolor="#f1f2f6", range=[0, max(plot_df["Deletion Frequency (%)"].max() * 1.15, 20)])
    )
    
    return fig

def plot_concordance_by_gene(gene_metrics_df):
    """
    Creates a horizontal bar chart showing accuracy/concordance rate per gene.
    """
    df = gene_metrics_df.copy()
    df["Concordance (%)"] = (df["Accuracy"] * 100).round(1)
    df = df.sort_values(by="Concordance (%)", ascending=True)
    
    fig = px.bar(
        df,
        y="Gene",
        x="Concordance (%)",
        orientation="h",
        color="Concordance (%)",
        color_continuous_scale="Tealgrn",
        title="Concordance Rate (%) by Gene",
        text="Concordance (%)"
    )
    
    fig.update_traces(textposition='inside', textfont_size=11, texttemplate='%{text}%')
    fig.update_layout(
        xaxis_title="Concordance Rate (%)",
        yaxis_title="Gene",
        coloraxis_showscale=False,
        margin=dict(l=40, r=40, t=60, b=40),
        height=450,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(gridcolor="#f1f2f6", range=[0, 105])
    )
    
    return fig

def plot_class_distribution(long_df, selected_genes=None):
    """
    Creates a stacked bar chart showing the breakdown of TP, TN, FP, FN
    for either all genes or specific selected genes.
    """
    df = long_df.copy()
    if selected_genes:
        df = df[df["Gene"].isin(selected_genes)]
        
    counts = df.groupby(["Gene", "Class"]).size().reset_index(name="Count")
    
    # Calculate percentage within each gene
    total_per_gene = df.groupby("Gene").size().reset_index(name="Total")
    counts = counts.merge(total_per_gene, on="Gene")
    counts["Percentage"] = (counts["Count"] / counts["Total"] * 100).round(1)
    
    # Map class labels
    counts["Class_Label"] = counts["Class"].map(COLOR_LABELS)
    
    # Order the classes for presentation
    class_order = ["TN", "TP", "FN", "FP"]
    
    fig = px.bar(
        counts,
        x="Gene",
        y="Percentage",
        color="Class",
        title="Comparison Breakdown (%) per Gene",
        color_discrete_map=COLOR_MAP,
        category_orders={"Class": class_order},
        hover_data=["Count", "Percentage"],
        labels={"Class": "Comparison Outcome"}
    )
    
    # Add human readable labels in the legend
    for i, trace in enumerate(fig.data):
        c_code = trace.name
        if c_code in COLOR_LABELS:
            trace.name = COLOR_LABELS[c_code]
            
    fig.update_layout(
        yaxis_title="Percentage (%)",
        xaxis_title="Gene",
        margin=dict(l=40, r=40, t=60, b=40),
        height=450,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        yaxis=dict(gridcolor="#f1f2f6", range=[0, 101]),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    
    return fig

def plot_cohort_heatmap(long_df, genes, sort_by="RegId"):
    """
    Creates an interactive clinical heatmap (grid) showing all patients (rows)
    and genes (columns), colored by TP, TN, FP, FN.
    This lets PIs quickly scan the cohort for systemic and case-specific patterns.
    """
    pivot_df = long_df.pivot(index="RegId", columns="Gene", values="Class")
    
    # Reindex columns to match the standard genes order
    pivot_df = pivot_df[genes]
    
    # Sort options
    if sort_by == "RegId":
        pivot_df = pivot_df.sort_index()
    elif sort_by == "Concordance":
        # Sort by concordance (count of TN + TP)
        concordance = long_df.groupby("RegId")["Concordant"].sum()
        pivot_df = pivot_df.loc[concordance.sort_values(ascending=False).index]
    elif sort_by == "Number of Deletions":
        # Sort by total deletions in DO
        deletions = long_df[long_df["DO_Status"] == "D"].groupby("RegId").size()
        # fill missing with 0
        all_patients = long_df["RegId"].unique()
        deletions = deletions.reindex(all_patients, fill_value=0)
        pivot_df = pivot_df.loc[deletions.sort_values(ascending=False).index]
        
    # Map classes to numeric values for heatmap plotting
    # TN=0, TP=1, FN=2, FP=3
    class_to_val = {"TN": 0, "TP": 1, "FN": 2, "FP": 3}
    val_to_class = {0: "TN", 1: "TP", 2: "FN", 3: "FP"}
    
    numeric_grid = pivot_df.replace(class_to_val).values
    y_labels = [f"Patient {rid}" for rid in pivot_df.index]
    x_labels = pivot_df.columns.tolist()
    
    # Define custom colorscale based on our COLOR_MAP
    # Values correspond to [0, 1, 2, 3] mapped to 0 to 1 normalized range
    colors = [COLOR_MAP["TN"], COLOR_MAP["TP"], COLOR_MAP["FN"], COLOR_MAP["FP"]]
    
    # Heatmap trace
    # To construct custom tooltips:
    hover_text = []
    for r_idx, rid in enumerate(pivot_df.index):
        row_text = []
        for c_idx, gene in enumerate(pivot_df.columns):
            outcome = pivot_df.loc[rid, gene]
            do_s = long_df[(long_df["RegId"] == rid) & (long_df["Gene"] == gene)]["DO_Status"].values[0]
            po_s = long_df[(long_df["RegId"] == rid) & (long_df["Gene"] == gene)]["PO_Status"].values[0]
            
            outcome_desc = COLOR_LABELS.get(outcome, outcome)
            
            txt = (
                f"<b>Patient ID:</b> {rid}<br>"
                f"<b>Gene:</b> {gene}<br>"
                f"<b>Comparison Outcome:</b> {outcome_desc}<br>"
                f"<b>Diagnostic (DO):</b> {'Deletion (D)' if do_s == 'D' else 'Neutral (N)'}<br>"
                f"<b>Pipeline (PO):</b> {'Deletion (D)' if po_s == 'D' else 'Neutral (N)'}"
            )
            row_text.append(txt)
        hover_text.append(row_text)
        
    # Heatmap figure
    fig = go.Figure(data=go.Heatmap(
        z=numeric_grid,
        x=x_labels,
        y=y_labels,
        colorscale=[
            [0.0, COLOR_MAP["TN"]], [0.25, COLOR_MAP["TN"]],
            [0.25, COLOR_MAP["TP"]], [0.5, COLOR_MAP["TP"]],
            [0.5, COLOR_MAP["FN"]], [0.75, COLOR_MAP["FN"]],
            [0.75, COLOR_MAP["FP"]], [1.0, COLOR_MAP["FP"]]
        ],
        showscale=False,
        text=hover_text,
        hoverinfo="text",
        xgap=2,
        ygap=2
    ))
    
    # Add dummy scatter traces to build a clean legend
    for cls_code, color in COLOR_MAP.items():
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode="markers",
            marker=dict(size=12, color=color, symbol="square"),
            name=COLOR_LABELS[cls_code],
            showlegend=True
        ))
        
    fig.update_layout(
        title=f"Clinical Cohort CNV Comparison Heatmap (Rows sorted by: {sort_by})",
        xaxis_title="Genes",
        yaxis_title="Patients",
        height=max(400, len(y_labels) * 15 + 100),  # Scale height with patient count
        margin=dict(l=100, r=40, t=60, b=40),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        yaxis=dict(autorange="reversed") # Show patient 1 at top
    )
    
    return fig

def plot_sensitivity_specificity(gene_metrics_df):
    """
    Creates a scatter plot showing Sensitivity vs Specificity for each gene.
    Marker size represents the F1-Score, and color represents Cohen's Kappa.
    """
    df = gene_metrics_df.copy()
    
    # Avoid division by zero/NaN issues for size
    df["F1_Size"] = df["F1_Score"].apply(lambda x: max(x, 0.1) * 30 + 10)
    
    fig = px.scatter(
        df,
        x="Sensitivity",
        y="Specificity",
        text="Gene",
        size="F1_Size",
        color="Kappa",
        color_continuous_scale="Viridis",
        title="Gene Validation Space: Sensitivity vs. Specificity (Marker size: F1-Score, Color: Kappa)",
        hover_name="Gene",
        hover_data={
            "Sensitivity": ":.3f",
            "Specificity": ":.3f",
            "Accuracy": ":.3f",
            "Precision": ":.3f",
            "F1_Score": ":.3f",
            "Kappa": ":.3f",
            "TP": True,
            "TN": True,
            "FP": True,
            "FN": True,
            "F1_Size": False
        }
    )
    
    # Adjust trace properties
    fig.update_traces(
        textposition="top center",
        marker=dict(line=dict(width=1, color="DarkSlateGrey"))
    )
    
    # Add reference diagonals or lines
    fig.add_shape(
        type="line", line=dict(dash="dash", color="lightgrey", width=1),
        x0=0.5, y0=0.5, x1=1, y1=1
    )
    
    fig.update_layout(
        xaxis=dict(title="Sensitivity (True Positive Rate)", range=[min(df["Sensitivity"].min() * 0.9, 0.5), 1.05], gridcolor="#f1f2f6"),
        yaxis=dict(title="Specificity (True Negative Rate)", range=[min(df["Specificity"].min() * 0.9, 0.5), 1.05], gridcolor="#f1f2f6"),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=40, r=40, t=60, b=40),
        height=500
    )
    
    return fig
