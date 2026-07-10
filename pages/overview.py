import streamlit as st
import pandas as pd
from utils.statistics import get_overall_metrics, get_gene_metrics, get_kappa_interpretation
from utils.plotting import plot_deletion_frequency, plot_concordance_by_gene

def run():
    # Page Header
    st.markdown("""
    <div class="portal-header">
        <h1>Clinical Overview Dashboard</h1>
        <p>High-level comparison of Pipeline Output (PO) against Diagnostic Output (DO) reference standard.</p>
    </div>
    """, unsafe_allow_html=True)

    long_df = st.session_state["long_df"]
    wide_df = st.session_state["wide_df"]
    genes = st.session_state["genes"]
    patients = st.session_state["patients"]

    # Calculate overall metrics
    metrics = get_overall_metrics(long_df)
    
    # Concordant patients count
    concordant_patients = wide_df[~wide_df["Has_Discrepancy"]].shape[0]
    pct_concordant_patients = (concordant_patients / len(patients)) * 100

    # Layout: Key Metrics
    col1, col2, col3, col4, col5 = st.columns(5)
    
    with col1:
        st.markdown(f"""
        <div class="metric-card">
            <h4 style="margin:0; color:#7f8c8d; font-size:0.9rem;">Cohort Size</h4>
            <p style="margin:5px 0 0 0; color:#2c3e50; font-size:1.8rem; font-weight:bold;">{len(patients)}</p>
            <span style="font-size:0.8rem; color:#7f8c8d;">Patients analyzed</span>
        </div>
        """, unsafe_allow_html=True)
        
    with col2:
        st.markdown(f"""
        <div class="metric-card">
            <h4 style="margin:0; color:#7f8c8d; font-size:0.9rem;">Overall Accuracy</h4>
            <p style="margin:5px 0 0 0; color:#10ac84; font-size:1.8rem; font-weight:bold;">{metrics['Accuracy']*100:.1f}%</p>
            <span style="font-size:0.8rem; color:#7f8c8d;">Gene-wise agreement</span>
        </div>
        """, unsafe_allow_html=True)
        
    with col3:
        st.markdown(f"""
        <div class="metric-card">
            <h4 style="margin:0; color:#7f8c8d; font-size:0.9rem;">Sensitivity (Recall)</h4>
            <p style="margin:5px 0 0 0; color:#2980b9; font-size:1.8rem; font-weight:bold;">{metrics['Sensitivity']*100:.1f}%</p>
            <span style="font-size:0.8rem; color:#7f8c8d;">True Deletion detection</span>
        </div>
        """, unsafe_allow_html=True)
        
    with col4:
        st.markdown(f"""
        <div class="metric-card">
            <h4 style="margin:0; color:#7f8c8d; font-size:0.9rem;">Specificity</h4>
            <p style="margin:5px 0 0 0; color:#8e44ad; font-size:1.8rem; font-weight:bold;">{metrics['Specificity']*100:.1f}%</p>
            <span style="font-size:0.8rem; color:#7f8c8d;">True Neutral detection</span>
        </div>
        """, unsafe_allow_html=True)
        
    with col5:
        st.markdown(f"""
        <div class="metric-card">
            <h4 style="margin:0; color:#7f8c8d; font-size:0.9rem;">Concordant Patients</h4>
            <p style="margin:5px 0 0 0; color:#27ae60; font-size:1.8rem; font-weight:bold;">{pct_concordant_patients:.1f}%</p>
            <span style="font-size:0.8rem; color:#27ae60; font-weight:500;">{concordant_patients}/{len(patients)} Patients</span>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # Main dashboard section: Summary stats & clinical context
    layout_col1, layout_col2 = st.columns([2, 3])
    
    with layout_col1:
        st.subheader("Clinical Summary")
        
        # Display Kappa
        kappa_val = metrics["Kappa"]
        kappa_interp = get_kappa_interpretation(kappa_val)
        
        st.info(f"""
        **Cohen's Kappa Agreement:** `{kappa_val:.3f}` ({kappa_interp})
        
        This metric measures the agreement between Diagnostic Output (DO) and Pipeline Output (PO), accounting for the agreement occurring by chance.
        """)
        
        # Discrepancy Breakdown description
        tp, tn, fp, fn = metrics["TP"], metrics["TN"], metrics["FP"], metrics["FN"]
        total_calls = tp + tn + fp + fn
        
        st.markdown(f"""
        ### Overall Call Breakdown
        Across all **{len(patients)} patients** and **{len(genes)} genes** (total of **{total_calls} calls**):
        - **True Negatives (TN):** `{tn}` ({tn/total_calls*100:.1f}%) - Both tests agree the gene is Neutral (N).
        - **True Positives (TP):** `{tp}` ({tp/total_calls*100:.1f}%) - Both tests agree the gene has a Deletion (D).
        - **False Positives (FP):** `{fp}` ({fp/total_calls*100:.1f}%) - The pipeline (PO) falsely detected a deletion.
        - **False Negatives (FN):** `{fn}` ({fn/total_calls*100:.1f}%) - The pipeline (PO) missed a deletion confirmed by DO.
        """)
        
        # Genes List
        st.markdown("### Genes Monitored")
        gene_badges = "".join([f'<span class="gene-badge">{g}</span>' for g in genes])
        st.markdown(f'<div style="line-height:2;">{gene_badges}</div>', unsafe_allow_html=True)

    with layout_col2:
        # Grouped bar chart DO vs PO
        fig_freq = plot_deletion_frequency(long_df)
        st.plotly_chart(fig_freq, width="stretch")

    st.markdown("<hr style='border: none; border-top: 1px solid #e1e8ed; margin: 30px 0;' />", unsafe_allow_html=True)

    # Detailed Gene Metrics Section
    st.subheader("Performance Metrics per Gene")
    gene_metrics_df = get_gene_metrics(long_df)
    
    # Format columns for display
    disp_metrics_df = gene_metrics_df.copy()
    disp_metrics_df["Accuracy"] = (disp_metrics_df["Accuracy"] * 100).round(1).astype(str) + "%"
    disp_metrics_df["Sensitivity"] = (disp_metrics_df["Sensitivity"] * 100).round(1).astype(str) + "%"
    disp_metrics_df["Specificity"] = (disp_metrics_df["Specificity"] * 100).round(1).astype(str) + "%"
    disp_metrics_df["Precision"] = (disp_metrics_df["Precision"] * 100).round(1).astype(str) + "%"
    disp_metrics_df["F1_Score"] = disp_metrics_df["F1_Score"].round(3)
    disp_metrics_df["Kappa"] = disp_metrics_df["Kappa"].round(3)
    
    # Apply style formatting
    st.dataframe(
        disp_metrics_df,
        column_config={
            "Gene": st.column_config.TextColumn("Gene Name", width="medium"),
            "TP": st.column_config.NumberColumn("True Pos (TP)"),
            "TN": st.column_config.NumberColumn("True Neg (TN)"),
            "FP": st.column_config.NumberColumn("False Pos (FP)"),
            "FN": st.column_config.NumberColumn("False Neg (FN)"),
            "Accuracy": st.column_config.TextColumn("Accuracy (Concordance)"),
            "Sensitivity": st.column_config.TextColumn("Sensitivity (Recall)"),
            "Specificity": st.column_config.TextColumn("Specificity"),
            "Precision": st.column_config.TextColumn("Precision (PPV)"),
            "F1_Score": st.column_config.NumberColumn("F1-Score"),
            "Kappa": st.column_config.NumberColumn("Cohen's Kappa")
        },
        width="stretch",
        hide_index=True
    )

    # Show a chart of concordance rate by gene
    fig_conc = plot_concordance_by_gene(gene_metrics_df)
    st.plotly_chart(fig_conc, width="stretch")

# Run if page file is executed directly (classic fallback)
if __name__ == "__main__":
    run()
elif __name__ == "pages.overview":
    run()
