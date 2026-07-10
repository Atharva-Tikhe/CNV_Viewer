import streamlit as st
import pandas as pd
import os
import sys

# Ensure utils directory is in path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from utils.loader import load_raw_data, get_long_comparison_df, get_patient_wide_comparison, get_genes, get_patients

# Set page config at the very beginning
st.set_page_config(
    page_title="Genomic CNV Concordance Portal",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling (Dark-themed glassmorphism elements, custom fonts, etc.)
st.markdown("""
<style>
    /* Global Styles */
    .stApp {
        background-color: #fafbfc;
    }
    
    /* Title and Header styling */
    .portal-header {
        background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
        padding: 24px;
        border-radius: 12px;
        color: white;
        margin-bottom: 24px;
        box-shadow: 0 4px 15px rgba(30, 60, 114, 0.15);
    }
    .portal-header h1 {
        margin: 0;
        font-family: 'Outfit', 'Inter', sans-serif;
        font-weight: 700;
        font-size: 2.2rem;
    }
    .portal-header p {
        margin: 8px 0 0 0;
        font-size: 1.05rem;
        opacity: 0.9;
    }
    
    /* Metrics Styling */
    .metric-card {
        background-color: white;
        border: 1px solid #e1e8ed;
        border-radius: 8px;
        padding: 16px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.02);
        transition: transform 0.2s, box-shadow 0.2s;
    }
    .metric-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 10px rgba(0,0,0,0.05);
    }
    
    /* Gene label styles */
    .gene-badge {
        display: inline-block;
        background-color: #ebf5fb;
        color: #2980b9;
        font-weight: 600;
        padding: 4px 8px;
        border-radius: 4px;
        font-size: 0.85rem;
        margin: 2px;
        border: 1px solid #d4e6f1;
    }
    
    /* Status Badge styling */
    .badge-tp { background-color: #e8f8f5; color: #117864; border: 1px solid #a3e4d7; font-weight: bold; padding: 2px 6px; border-radius: 4px; }
    .badge-tn { background-color: #f8f9f9; color: #5d6d7e; border: 1px solid #d5dbdb; font-weight: bold; padding: 2px 6px; border-radius: 4px; }
    .badge-fp { background-color: #fdebd0; color: #d35400; border: 1px solid #f5cba7; font-weight: bold; padding: 2px 6px; border-radius: 4px; }
    .badge-fn { background-color: #fadbd8; color: #78281f; border: 1px solid #f1948a; font-weight: bold; padding: 2px 6px; border-radius: 4px; }
</style>
""", unsafe_allow_html=True)

# Initialize Session State
if "raw_df" not in st.session_state:
    try:
        raw_df = load_raw_data()
        st.session_state["raw_df"] = raw_df
        st.session_state["long_df"] = get_long_comparison_df(raw_df)
        st.session_state["wide_df"] = get_patient_wide_comparison(raw_df)
        st.session_state["genes"] = get_genes(raw_df)
        st.session_state["patients"] = get_patients(raw_df)
    except Exception as e:
        st.error(f"Error loading data: {e}")

# Notes database in session state (so PIs can annotate cases)
if "patient_notes" not in st.session_state:
    st.session_state["patient_notes"] = {}

# Active filter settings stored globally
if "active_filters" not in st.session_state:
    st.session_state["active_filters"] = {
        "genes": st.session_state.get("genes", []),
        "patients": st.session_state.get("patients", []),
        "concordance": "All",
        "gene_deletions": {}
    }

# Multipage Navigation Setup
# Standard folder layout pages are:
# pages/overview.py, pages/patient_browser.py, pages/vis.py, pages/stats.py, pages/filters.py, pages/export.py
pages = {
    "Dashboard": [
        st.Page("pages/overview.py", title="Overview Dashboard", icon="📊"),
        st.Page("pages/patient_browser.py", title="Patient Case Browser", icon="👤"),
    ],
    "Analytics": [
        st.Page("pages/vis.py", title="Interactive Charts", icon="📈"),
        st.Page("pages/stats.py", title="Validation Statistics", icon="🔬"),
    ],
    "Explorer": [
        st.Page("pages/filters.py", title="Deletion Pattern Explorer", icon="🔍"),
        st.Page("pages/export.py", title="Export & Reporting", icon="💾"),
    ]
}

# Add general metadata in sidebar
st.sidebar.markdown("""
<div style="text-align: center; padding: 10px 0;">
    <h2 style="margin: 0; color: #1e3c72; font-size: 1.4rem;">CNV Explorer</h2>
    <p style="margin: 0; font-size: 0.85rem; color: #7f8c8d;">Version 1.0 (Illumina | Batch 6)</p>
</div>
<hr style="margin: 10px 0; border: none; border-top: 1px solid #e1e8ed;" />
""", unsafe_allow_html=True)

# Navigation execution
if hasattr(st, "navigation"):
    pg = st.navigation(pages)
    pg.run()
else:
    # Classic sidebar fallback for older Streamlit versions
    st.sidebar.title("Navigation")
    page_names = ["Overview Dashboard", "Patient Case Browser", "Interactive Charts", "Validation Statistics", "Deletion Pattern Explorer", "Export & Reporting"]
    selected_page = st.sidebar.radio("Go to", page_names)
    
    # Fallback import & run
    if selected_page == "Overview Dashboard":
        import pages.overview as page_module
    elif selected_page == "Patient Case Browser":
        import pages.patient_browser as page_module
    elif selected_page == "Interactive Charts":
        import pages.vis as page_module
    elif selected_page == "Validation Statistics":
        import pages.stats as page_module
    elif selected_page == "Deletion Pattern Explorer":
        import pages.filters as page_module
    elif selected_page == "Export & Reporting":
        import pages.export as page_module
