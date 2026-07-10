import streamlit as st
import pandas as pd
import json
import os

# Path for saving notes persistently
NOTES_FILE = "data/validation_notes.json"

def load_persistent_notes():
    """Loads notes from disk if they exist, else returns empty dict."""
    if os.path.exists(NOTES_FILE):
        try:
            with open(NOTES_FILE, "r") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}

def save_persistent_notes(notes):
    """Saves notes to disk."""
    try:
        os.makedirs(os.path.dirname(NOTES_FILE), exist_ok=True)
        with open(NOTES_FILE, "w") as f:
            json.dump(notes, f, indent=4)
        return True
    except Exception as e:
        st.error(f"Error saving notes to file: {e}")
        return False

def run():
    # Header
    st.markdown("""
    <div class="portal-header">
        <h1>👤 Patient Case Browser</h1>
        <p>Explore gene-by-gene results for individual patients, document review comments, and log clinical validations.</p>
    </div>
    """, unsafe_allow_html=True)

    long_df = st.session_state["long_df"]
    wide_df = st.session_state["wide_df"]
    genes = st.session_state["genes"]
    patients = st.session_state["patients"]

    # Load persistent notes
    if "persistent_notes" not in st.session_state:
        st.session_state["persistent_notes"] = load_persistent_notes()
        
    persistent_notes = st.session_state["persistent_notes"]

    # Sidebar / Selection panel
    st.sidebar.subheader("Patient Selection")
    
    # Quick filter in sidebar for patients with discrepancies
    only_discrepant = st.sidebar.checkbox("Show only patients with discrepancies", value=False)
    
    if only_discrepant:
        filtered_patients = wide_df[wide_df["Has_Discrepancy"]]["RegId"].tolist()
        if not filtered_patients:
            st.sidebar.warning("No patients with discrepancies found.")
            filtered_patients = patients
    else:
        filtered_patients = patients

    selected_rid = st.sidebar.selectbox("Select Patient (RegID)", filtered_patients)

    if selected_rid:
        # Extract patient data
        p_long = long_df[long_df["RegId"] == selected_rid].copy()
        p_wide = wide_df[wide_df["RegId"] == selected_rid].iloc[0]
        
        # Patient overview summary
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h4 style="margin:0; color:#7f8c8d; font-size:0.9rem;">Patient ID</h4>
                <p style="margin:5px 0 0 0; color:#2c3e50; font-size:1.6rem; font-weight:bold;">{selected_rid}</p>
                <span style="font-size:0.8rem; color:#7f8c8d;">Registration Identifier</span>
            </div>
            """, unsafe_allow_html=True)
            
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h4 style="margin:0; color:#7f8c8d; font-size:0.9rem;">Concordance Rate</h4>
                <p style="margin:5px 0 0 0; color:#10ac84; font-size:1.6rem; font-weight:bold;">{p_wide['Concordance_Rate']:.1f}%</p>
                <span style="font-size:0.8rem; color:#7f8c8d;">{int(p_wide['Concordant_Genes'])} / {int(p_wide['Total_Genes'])} genes concordant</span>
            </div>
            """, unsafe_allow_html=True)
            
        with col3:
            st.markdown(f"""
            <div class="metric-card">
                <h4 style="margin:0; color:#7f8c8d; font-size:0.9rem;">Deletions (DO vs PO)</h4>
                <p style="margin:5px 0 0 0; color:#2980b9; font-size:1.6rem; font-weight:bold;">{int(p_wide['TP'] + p_wide['FN'])} vs {int(p_wide['TP'] + p_wide['FP'])}</p>
                <span style="font-size:0.8rem; color:#7f8c8d;">Diagnostic vs Pipeline</span>
            </div>
            """, unsafe_allow_html=True)
            
        with col4:
            discrepancies = int(p_wide['FP'] + p_wide['FN'])
            status_color = "#ee5253" if discrepancies > 0 else "#10ac84"
            status_text = f"{discrepancies} Discrepancy" if discrepancies == 1 else f"{discrepancies} Discrepancies"
            
            st.markdown(f"""
            <div class="metric-card">
                <h4 style="margin:0; color:#7f8c8d; font-size:0.9rem;">Discrepancy Status</h4>
                <p style="margin:5px 0 0 0; color:{status_color}; font-size:1.6rem; font-weight:bold;">{status_text}</p>
                <span style="font-size:0.8rem; color:#7f8c8d;">{"Needs clinical review" if discrepancies > 0 else "Fully concordant"}</span>
            </div>
            """, unsafe_allow_html=True)

        st.markdown("<br>", unsafe_allow_html=True)

        # Content split: Left is CNV comparison, Right is review notes
        c_left, c_right = st.columns([5, 4])
        
        with c_left:
            st.subheader("Gene-by-Gene CNV Comparison")
            
            # Build clean display table
            disp_df = p_long[["Gene", "DO_Status", "PO_Status", "Class", "Class_Label"]].copy()
            disp_df.columns = ["Gene Name", "Diagnostic (DO)", "Pipeline (PO)", "Class", "Review Outcome"]
            
            # Map status values to descriptive text
            status_map = {"D": "Deletion (D)", "N": "Neutral (N)"}
            disp_df["Diagnostic (DO)"] = disp_df["Diagnostic (DO)"].map(status_map)
            disp_df["Pipeline (PO)"] = disp_df["Pipeline (PO)"].map(status_map)
            
            # Styling function
            def style_cells(data):
                styles = pd.DataFrame('', index=data.index, columns=data.columns)
                for idx, row in data.iterrows():
                    cls = row['Class']
                    if cls == 'TP': # True Positive (Match Deletion)
                        color = '#e8f8f5'
                        text_color = '#117864'
                    elif cls == 'TN': # True Negative (Match Neutral)
                        color = '#f8f9f9'
                        text_color = '#5d6d7e'
                    elif cls == 'FP': # False Positive
                        color = '#fdebd0'
                        text_color = '#d35400'
                    elif cls == 'FN': # False Negative
                        color = '#fadbd8'
                        text_color = '#78281f'
                    else:
                        color = ''
                        text_color = ''
                        
                    styles.loc[idx, 'Review Outcome'] = f'background-color: {color}; color: {text_color}; font-weight: bold;'
                    styles.loc[idx, 'Gene Name'] = 'font-weight: bold;'
                    
                    # Style DO
                    if 'Deletion' in row['Diagnostic (DO)']:
                        styles.loc[idx, 'Diagnostic (DO)'] = 'color: #e74c3c; font-weight: bold;'
                    else:
                        styles.loc[idx, 'Diagnostic (DO)'] = 'color: #34495e;'
                        
                    # Style PO
                    if 'Deletion' in row['Pipeline (PO)']:
                        styles.loc[idx, 'Pipeline (PO)'] = 'color: #e74c3c; font-weight: bold;'
                    else:
                        styles.loc[idx, 'Pipeline (PO)'] = 'color: #34495e;'
                        
                return styles
            
            styled_df = disp_df.style.apply(style_cells, axis=None)
            
            st.dataframe(
                styled_df,
                column_config={
                    "Class": None, # Hides the helper class column
                },
                width="stretch",
                hide_index=True
            )
            
        with c_right:
            st.subheader("PI Clinical Validation Notes")
            
            # Check if notes exist for this patient
            p_note_key = str(selected_rid)
            existing_note = persistent_notes.get(p_note_key, {})
            
            e_status = existing_note.get("status", "Pending Review")
            e_comment = existing_note.get("comment", "")
            e_reviewer = existing_note.get("reviewer", "")
            
            # Notes Form
            with st.form("validation_notes_form"):
                review_status = st.selectbox(
                    "Validation Status",
                    [
                        "Pending Review", 
                        "Concordant - Approved", 
                        "Discrepancy - PO Falsely Called (FP)", 
                        "Discrepancy - PO Missed (FN)", 
                        "Discrepancy - Reference DO Doubtful",
                        "Technical Failure / Resequencing Needed"
                    ],
                    index=[
                        "Pending Review", 
                        "Concordant - Approved", 
                        "Discrepancy - PO Falsely Called (FP)", 
                        "Discrepancy - PO Missed (FN)", 
                        "Discrepancy - Reference DO Doubtful",
                        "Technical Failure / Resequencing Needed"
                    ].index(e_status) if e_status in [
                        "Pending Review", 
                        "Concordant - Approved", 
                        "Discrepancy - PO Falsely Called (FP)", 
                        "Discrepancy - PO Missed (FN)", 
                        "Discrepancy - Reference DO Doubtful",
                        "Technical Failure / Resequencing Needed"
                    ] else 0
                )
                
                reviewer_name = st.text_input("Reviewer Initials / Name", value=e_reviewer, placeholder="e.g. Dr. Jane Doe")
                review_comment = st.text_area("Clinical Notes & Remarks", value=e_comment, placeholder="Enter any validation remarks, data quality comments, or reasoning here...", height=120)
                
                submit_button = st.form_submit_button("Save Review Notes", width="stretch")
                
                if submit_button:
                    # Update notes dictionary
                    persistent_notes[p_note_key] = {
                        "status": review_status,
                        "comment": review_comment,
                        "reviewer": reviewer_name,
                        "last_updated": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
                    }
                    st.session_state["persistent_notes"] = persistent_notes
                    
                    # Save to file
                    if save_persistent_notes(persistent_notes):
                        st.success("Review notes successfully saved!")
                        st.rerun()

            # Display currently logged note if it exists
            if p_note_key in persistent_notes:
                st.markdown("### Logged Review Detail")
                note_detail = persistent_notes[p_note_key]
                st.markdown(f"""
                * **Status:** `{note_detail['status']}`
                * **Reviewer:** `{note_detail['reviewer'] or 'Not specified'}`
                * **Last Updated:** `{note_detail['last_updated']}`
                * **Comment:**
                > {note_detail['comment'] or '*No comment entered.*'}
                """)
                
# Run if page file is executed directly (classic fallback)
if __name__ == "__main__":
    run()
elif __name__ == "pages.patient_browser":
    run()
