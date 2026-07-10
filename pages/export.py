import streamlit as st
import pandas as pd
import json
import io
import os

NOTES_FILE = "data/validation_notes.json"

def get_reviewer_summary(notes):
    """Summarizes reviewer statuses."""
    if not notes:
        return {}
    statuses = {}
    for rid, note in notes.items():
        st_val = note.get("status", "Pending Review")
        statuses[st_val] = statuses.get(st_val, 0) + 1
    return statuses

def run():
    st.markdown("""
    <div class="portal-header">
        <h1>💾 Export & Reporting</h1>
        <p>Generate clinical reports, download filtered datasets, and export PI validation review logs.</p>
    </div>
    """, unsafe_allow_html=True)

    long_df = st.session_state["long_df"]
    wide_df = st.session_state["wide_df"]
    raw_df = st.session_state["raw_df"]
    genes = st.session_state["genes"]
    patients = st.session_state["patients"]

    # Load active filters from session state to filter the downloaded data
    active_filters = st.session_state.get("active_filters", {})
    selected_genes = active_filters.get("genes", genes)
    concordance = active_filters.get("concordance", "All")
    
    # Reload notes
    if os.path.exists(NOTES_FILE):
        try:
            with open(NOTES_FILE, "r") as f:
                notes = json.load(f)
        except Exception:
            notes = {}
    else:
        notes = {}

    st.subheader("Data Export Center")
    st.markdown("Download datasets in Excel or CSV formats for clinical reporting or downstream analysis.")

    # Prepare datasets
    filtered_wide = wide_df
    
    # 1. Excel Exporter
    # We will build an Excel file in-memory using openpyxl and offer it as a download
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        raw_df.to_excel(writer, sheet_name='Raw Data', index=False)
        long_df.to_excel(writer, sheet_name='Long Comparisons', index=False)
        wide_df.to_excel(writer, sheet_name='Patient Summaries', index=False)
        
        # Include notes sheet if they exist
        if notes:
            notes_list = []
            for rid, note in notes.items():
                notes_list.append({
                    "RegId": rid,
                    "Status": note.get("status"),
                    "Reviewer": note.get("reviewer"),
                    "Comment": note.get("comment"),
                    "Last Updated": note.get("last_updated")
                })
            pd.DataFrame(notes_list).to_excel(writer, sheet_name='PI Notes Logs', index=False)
            
    excel_data = buffer.getvalue()

    # Layout for downloads
    c_xlsx, c_csv_long, c_csv_wide = st.columns(3)
    
    with c_xlsx:
        st.markdown("### Combined Excel Report")
        st.caption("Contains raw, long-comparison, patient-wide summaries, and review sheets in a single file.")
        st.download_button(
            label="📥 Download Excel Workbook",
            data=excel_data,
            file_name="CNV_Concordance_Report.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            width="stretch"
        )

    with c_csv_long:
        st.markdown("### Long Comparison CSV")
        st.caption("Detailed call-by-call classification (TP, TN, FP, FN) for every gene and patient.")
        csv_long_data = long_df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Long CSV",
            data=csv_long_data,
            file_name="cnv_long_comparisons.csv",
            mime="text/csv",
            width="stretch"
        )

    with c_csv_wide:
        st.markdown("### Patient Summary CSV")
        st.caption("One row per patient including individual gene comparisons and overall patient concordance metrics.")
        csv_wide_data = wide_df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Wide CSV",
            data=csv_wide_data,
            file_name="cnv_patient_summaries.csv",
            mime="text/csv",
            width="stretch"
        )

    st.markdown("<hr style='border: none; border-top: 1px solid #e1e8ed; margin: 30px 0;' />", unsafe_allow_html=True)

    # Summary Report Builder for PI
    st.subheader("Clinical Validation Report")
    
    # Generate automated report
    from utils.statistics import get_overall_metrics, get_gene_metrics
    metrics = get_overall_metrics(long_df)
    
    status_summary = get_reviewer_summary(notes)
    pending_count = len(patients) - len(notes)
    
    report_text = f"""# Clinical CNV Validation Report
Generated on: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}
Dataset: Illumina | Batch 6

## Cohort Summary
* **Total Cohort Size:** {len(patients)} Patients
* **Number of Genes Analyzed:** {len(genes)}
* **Overall Concordance Rate:** {metrics['Accuracy']*100:.2f}%
* **Overall Deletion Sensitivity (Recall):** {metrics['Sensitivity']*100:.2f}%
* **Overall Deletion Specificity:** {metrics['Specificity']*100:.2f}%
* **Cohen's Kappa Score:** {metrics['Kappa']:.3f}

## PI Validation Review Status
* **Pending Review:** {pending_count} cases
* **Approved Concordant:** {status_summary.get('Concordant - Approved', 0)} cases
* **Confirmed Pipeline False Positives (FP):** {status_summary.get('Discrepancy - PO Falsely Called (FP)', 0)} cases
* **Confirmed Pipeline False Negatives (FN):** {status_summary.get('Discrepancy - PO Missed (FN)', 0)} cases
* **Doubtful Reference (DO):** {status_summary.get('Discrepancy - Reference DO Doubtful', 0)} cases
* **Technical Failure / Re-run required:** {status_summary.get('Technical Failure / Resequencing Needed', 0)} cases

## Top Discrepant Genes (Ordered by lowest Kappa)
"""
    # Find discrepant genes
    gene_metrics = get_gene_metrics(long_df)
    worst_genes = gene_metrics.sort_values(by="Kappa").head(3)
    
    for idx, row in worst_genes.iterrows():
        report_text += f"* **{row['Gene']}:** Kappa = `{row['Kappa']:.3f}`, Sensitivity = `{row['Sensitivity']*100:.1f}%`, Specificity = `{row['Specificity']*100:.1f}%` (FP: `{row['FP']}`, FN: `{row['FN']}`)\n"

    report_text += """
## Notes & Remarks
Report generated automatically by CNV Validation Portal. Use this report for clinical regulatory documentation and validation audits.
"""

    st.markdown(report_text)
    
    st.download_button(
        label="📥 Download Text Report (.md)",
        data=report_text.encode('utf-8'),
        file_name="CNV_Clinical_Report.md",
        mime="text/markdown",
        width="content"
    )

# Run if page file is executed directly (classic fallback)
if __name__ == "__main__":
    run()
elif __name__ == "pages.export":
    run()
