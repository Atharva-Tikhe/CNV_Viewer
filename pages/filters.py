import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from utils.loader import get_genes

def compute_deletion_patterns(raw_df, source="DO"):
    """
    Computes deletion frequencies and conditional probabilities (co-occurrence)
    for each gene in the cohort under a specific output type (DO or PO).
    """
    genes = get_genes(raw_df)
    
    # Filter for the chosen source (DO or PO)
    df_source = raw_df[raw_df["Source"] == source].copy()
    
    # Convert N -> 0, D -> 1
    binary_df = df_source[genes].copy()
    for g in genes:
        binary_df[g] = binary_df[g].map({"D": 1, "N": 0}).fillna(0)
        
    total_samples = len(df_source)
    
    # 1. Calculate deletion counts and percentages
    del_counts = binary_df.sum()
    del_pct = (del_counts / total_samples * 100).round(1)
    
    summary_df = pd.DataFrame({
        "Gene": genes,
        "Deleted Samples": del_counts.values,
        "Deletion Rate (%)": del_pct.values
    }).sort_values(by="Deletion Rate (%)", ascending=False)
    
    # 2. Compute co-occurrence matrix (both deleted)
    # Dot product of transpose with itself gives intersection counts
    intersection = binary_df.T.dot(binary_df)
    
    # 3. Compute conditional probability: P(Col Gene is deleted | Row Gene is deleted)
    # P(Col | Row) = (Row and Col) / Row
    # Divide each row of intersection by the deletion count of that row gene
    cond_prob = intersection.div(del_counts, axis=0) * 100
    cond_prob = cond_prob.fillna(0).round(1)
    
    return summary_df, intersection, cond_prob, total_samples

def run():
    st.markdown("""
    <div class="portal-header">
        <h1>Deletion Pattern Explorer</h1>
        <p>Analyze gene deletion frequencies and discover co-occurrence patterns across samples.</p>
    </div>
    """, unsafe_allow_html=True)

    raw_df = st.session_state["raw_df"]
    genes = st.session_state["genes"]

    # Source Selection
    st.markdown("### Select Data Source for Analysis")
    source_choice = st.radio(
        "Analyze patterns in:",
        ["Diagnostic Output (DO - Reference Standard)", "Pipeline Output (PO - Pipeline Calls)"],
        horizontal=True,
        help="Analyze deletion patterns either in the diagnostic standard or in the pipeline's output."
    )
    
    source = "DO" if "Diagnostic" in source_choice else "PO"

    # Compute patterns
    summary_df, intersection, cond_prob, total_samples = compute_deletion_patterns(raw_df, source=source)

    # Display Layout
    st.subheader("Gene Deletion Frequencies")
    st.markdown(f"Total samples analyzed: **{total_samples}**")
    
    # Create side-by-side display: bar chart and table
    col_chart, col_table = st.columns([3, 2])
    
    with col_chart:
        # Plot deletion frequencies
        fig_freq = px.bar(
            summary_df,
            x="Gene",
            y="Deletion Rate (%)",
            color="Deletion Rate (%)",
            color_continuous_scale="Reds",
            title=f"Gene Deletion Rates in {source}",
            text="Deletion Rate (%)"
        )
        fig_freq.update_traces(textposition="outside", texttemplate="%{text}%")
        fig_freq.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            yaxis=dict(gridcolor="#f1f2f6", range=[0, 115]),
            coloraxis_showscale=False,
            height=380,
            margin=dict(l=40, r=40, t=50, b=40)
        )
        st.plotly_chart(fig_freq, width="stretch")
        
    with col_table:
        st.markdown("<br>", unsafe_allow_html=True)
        st.dataframe(
            summary_df,
            column_config={
                "Gene": st.column_config.TextColumn("Gene Name"),
                "Deleted Samples": st.column_config.NumberColumn("Deleted Count"),
                "Deletion Rate (%)": st.column_config.NumberColumn("Rate (%)", format="%.1f%%")
            },
            width="stretch",
            hide_index=True
        )

    st.markdown("<hr style='border: none; border-top: 1px solid #e1e8ed; margin: 30px 0;' />", unsafe_allow_html=True)

    # Co-occurrence Section
    st.subheader("Deletion Co-occurrence & Conditional Patterns")
    st.markdown("""
    Explore dependencies between gene deletions. Select a primary gene to see the probability of other genes being deleted in the same sample.
    """)

    # Filter out genes that have 0 deletions to avoid divide by zero text issues
    active_genes = summary_df[summary_df["Deleted Samples"] > 0]["Gene"].tolist()
    
    if not active_genes:
        st.warning("No deletions detected in the dataset for the selected source.")
        return

    col_select, col_info = st.columns([1, 2])
    with col_select:
        primary_gene = st.selectbox(
            "Select Primary Gene",
            active_genes,
            help="Select the gene that you assume is deleted to study co-deletions."
        )
        
    primary_del_count = int(summary_df[summary_df["Gene"] == primary_gene]["Deleted Samples"].values[0])
    
    with col_info:
        st.markdown(f"""
        <div style="background-color: #ebf5fb; padding: 12px 18px; border-radius: 8px; border-left: 5px solid #2980b9; margin-top: 5px;">
            <p style="margin: 0; font-size: 0.95rem; color: #2c3e50;">
                Analyzing samples where <b>{primary_gene}</b> is deleted (<b>{primary_del_count} samples</b>).
            </p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # Calculate conditional probabilities for primary gene
    # Row corresponding to primary_gene in cond_prob dataframe
    gene_conds = cond_prob.loc[primary_gene].reset_index()
    gene_conds.columns = ["Gene", "Conditional Deletion Probability (%)"]
    
    # Filter out the primary gene itself (which is always 100%)
    gene_conds_other = gene_conds[gene_conds["Gene"] != primary_gene].sort_values(
        by="Conditional Deletion Probability (%)", ascending=False
    )
    
    # Add count columns for clarity
    # Both deleted count
    both_counts = []
    for g in gene_conds_other["Gene"]:
        both_counts.append(int(intersection.loc[primary_gene, g]))
    gene_conds_other["Co-deleted Count"] = both_counts
    
    # Layout for conditional probability details
    c_prob_chart, c_text = st.columns([3, 2])
    
    with c_prob_chart:
        fig_cond = px.bar(
            gene_conds_other,
            x="Gene",
            y="Conditional Deletion Probability (%)",
            title=f"Probability of Co-deletion Given {primary_gene} is Deleted",
            color="Conditional Deletion Probability (%)",
            color_continuous_scale="Oranges",
            text="Conditional Deletion Probability (%)"
        )
        fig_cond.update_traces(textposition="outside", texttemplate="%{text}%")
        fig_cond.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            yaxis=dict(gridcolor="#f1f2f6", range=[0, 115]),
            coloraxis_showscale=False,
            height=350,
            margin=dict(l=40, r=40, t=50, b=40)
        )
        st.plotly_chart(fig_cond, width="stretch")
        
    with c_text:
        st.markdown(f"#### Co-occurrence Details for {primary_gene}")
        for _, r in gene_conds_other.iterrows():
            other_g = r["Gene"]
            prob = r["Conditional Deletion Probability (%)"]
            co_cnt = int(r["Co-deleted Count"])
            
            if prob > 0:
                st.markdown(f"* When **{primary_gene}** is deleted, **{other_g}** is also deleted **{prob:.1f}%** of the time ({co_cnt}/{primary_del_count} samples).")
            else:
                st.markdown(f"* **{other_g}** is never co-deleted with **{primary_gene}**.")

    st.markdown("<hr style='border: none; border-top: 1px solid #e1e8ed; margin: 30px 0;' />", unsafe_allow_html=True)

    # 3. Heatmap of Conditional Deletion Probability
    st.subheader("Cohort Conditional Deletion Probability Matrix")
    st.markdown("""
    This matrix shows the conditional probability $P(\text{Column Gene deleted} \mid \text{Row Gene deleted})$.
    Hover over each cell to read the exact co-deletion probability and sample overlaps.
    """)
    
    # Re-order matrix to match summary_df gene order (highest deletion rate first) for readability
    sorted_genes = summary_df["Gene"].tolist()
    sorted_cond_prob = cond_prob.loc[sorted_genes, sorted_genes]
    
    # Custom hover text
    hover_t = []
    for r_gene in sorted_genes:
        row_text = []
        r_del_cnt = int(summary_df[summary_df["Gene"] == r_gene]["Deleted Samples"].values[0])
        for c_gene in sorted_genes:
            p_val = sorted_cond_prob.loc[r_gene, c_gene]
            both_cnt = int(intersection.loc[r_gene, c_gene])
            
            if r_gene == c_gene:
                txt = f"<b>{r_gene}</b><br>Total Deleted: {r_del_cnt} samples"
            else:
                txt = (
                    f"<b>Given:</b> {r_gene} is deleted ({r_del_cnt} samples)<br>"
                    f"<b>Then:</b> {c_gene} is also deleted ({both_cnt} samples)<br>"
                    f"<b>Probability:</b> {p_val:.1f}%"
                )
            row_text.append(txt)
        hover_t.append(row_text)
        
    fig_mat = go.Figure(data=go.Heatmap(
        z=sorted_cond_prob.values,
        x=sorted_genes,
        y=sorted_genes,
        colorscale="YlOrRd",
        text=hover_t,
        hoverinfo="text",
        xgap=1,
        ygap=1
    ))
    
    fig_mat.update_layout(
        xaxis_title="Co-deleted Gene (Column)",
        yaxis_title="Given Deleted Gene (Row)",
        height=500,
        margin=dict(l=80, r=40, t=50, b=40),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        yaxis=dict(autorange="reversed")
    )
    
    st.plotly_chart(fig_mat, width="stretch")

# Run if page file is executed directly (classic fallback)
if __name__ == "__main__":
    run()
elif __name__ == "pages.filters":
    run()
