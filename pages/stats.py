import streamlit as st
import pandas as pd
import numpy as np
from utils.statistics import get_overall_metrics, get_gene_metrics, get_kappa_interpretation, calculate_metrics

def get_mcnemar_p_value(fp, fn):
    """
    Calculates the exact binomial p-value (or corrected chi-square p-value) for McNemar's test.
    If fp + fn < 25, uses the exact binomial distribution (2-tailed).
    Else uses the chi-squared distribution with 1 degree of freedom (with continuity correction).
    """
    total_discordant = fp + fn
    if total_discordant == 0:
        return 1.0, "No discrepancies"
        
    # Exact binomial p-value
    from scipy.stats import binom
    # Under H0, the discordant pairs are split 50/50.
    # So we compute the probability of getting min(fp, fn) or fewer successes in total_discordant trials with p=0.5
    # and multiply by 2 for a two-tailed test.
    exact_p = 2 * binom.cdf(min(fp, fn), total_discordant, 0.5)
    exact_p = min(exact_p, 1.0)
    
    # Asymptotic Chi-Square with continuity correction
    chi2_stat = (abs(fp - fn) - 1)**2 / total_discordant if total_discordant > 0 else 0
    
    # Return both
    if total_discordant < 25:
        p_val = exact_p
        method = "McNemar's Exact Binomial Test (recommended for small samples)"
    else:
        # Using chi2 approximation. We can calculate or just stick to exact binomial which is always valid and accurate!
        # Scipy is installed in venv, so we can use scipy.stats.chi2
        from scipy.stats import chi2
        p_val = chi2.sf(chi2_stat, 1)
        method = "McNemar's Chi-squared Test (with continuity correction)"
        
    return p_val, method

def run():
    st.markdown("""
    <div class="portal-header">
        <h1>🔬 Validation Statistics</h1>
        <p>In-depth statistical validation, confusion matrices, and agreement tests (Cohen's Kappa & McNemar's Test).</p>
    </div>
    """, unsafe_allow_html=True)

    long_df = st.session_state["long_df"]
    genes = st.session_state["genes"]
    
    # Selector for Analysis Level
    st.markdown("### Analysis Level")
    level = st.selectbox("Select Target for Statistical Review", ["Overall Cohort (All Genes Combined)"] + genes)
    
    if level == "Overall Cohort (All Genes Combined)":
        metrics = get_overall_metrics(long_df)
        title_suffix = "Overall Cohort"
    else:
        gene_data = long_df[long_df["Gene"] == level]
        tp = (gene_data["Class"] == "TP").sum()
        tn = (gene_data["Class"] == "TN").sum()
        fp = (gene_data["Class"] == "FP").sum()
        fn = (gene_data["Class"] == "FN").sum()
        metrics = calculate_metrics(tp, tn, fp, fn)
        title_suffix = f"Gene: {level}"

    tp, tn, fp, fn = metrics["TP"], metrics["TN"], metrics["FP"], metrics["FN"]

    # Show layout: Confusion Matrix on left, Detailed Metrics on right
    col_cm, col_metrics = st.columns([1, 1.2])
    
    with col_cm:
        st.subheader("Confusion Matrix")
        
        # Display Confusion Matrix as a beautiful HTML Table
        st.markdown(f"""
        <div style="background-color: white; border: 1px solid #e1e8ed; border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.01);">
            <table style="width:100%; text-align: center; border-collapse: collapse; font-family: sans-serif;">
                <tr>
                    <td colspan="2" rowspan="2" style="border: none;"></td>
                    <td colspan="2" style="background-color: #3498db; color: white; font-weight: bold; padding: 10px; border-radius: 6px 6px 0 0;">Pipeline Output (PO)</td>
                </tr>
                <tr>
                    <td style="background-color: #ebf5fb; color: #2980b9; font-weight: bold; padding: 10px; width: 40%; border: 1px solid #d4e6f1;">Deletion (D)</td>
                    <td style="background-color: #ebf5fb; color: #2980b9; font-weight: bold; padding: 10px; width: 40%; border: 1px solid #d4e6f1;">Neutral (N)</td>
                </tr>
                <tr>
                    <td rowspan="2" style="background-color: #34495e; color: white; font-weight: bold; padding: 10px; writing-mode: vertical-rl; transform: rotate(180deg); border-radius: 6px 0 0 6px;">Diagnostic (DO)</td>
                    <td style="background-color: #f2f4f4; color: #34495e; font-weight: bold; padding: 10px; border: 1px solid #d5dbdb;">Deletion (D)</td>
                    <td style="background-color: #e8f8f5; color: #117864; font-weight: bold; padding: 15px; border: 2px solid #a3e4d7; font-size: 1.2rem;">
                        {tp}<br><span style="font-size: 0.8rem; font-weight: normal; color: #16a085;">True Positive (TP)</span>
                    </td>
                    <td style="background-color: #fadbd8; color: #78281f; font-weight: bold; padding: 15px; border: 2px solid #f1948a; font-size: 1.2rem;">
                        {fn}<br><span style="font-size: 0.8rem; font-weight: normal; color: #c0392b;">False Negative (FN)</span>
                    </td>
                </tr>
                <tr>
                    <td style="background-color: #f2f4f4; color: #34495e; font-weight: bold; padding: 10px; border: 1px solid #d5dbdb;">Neutral (N)</td>
                    <td style="background-color: #fdebd0; color: #d35400; font-weight: bold; padding: 15px; border: 2px solid #f5cba7; font-size: 1.2rem;">
                        {fp}<br><span style="font-size: 0.8rem; font-weight: normal; color: #d35400;">False Positive (FP)</span>
                    </td>
                    <td style="background-color: #f8f9f9; color: #5d6d7e; font-weight: bold; padding: 15px; border: 2px solid #d5dbdb; font-size: 1.2rem;">
                        {tn}<br><span style="font-size: 0.8rem; font-weight: normal; color: #7f8c8d;">True Negative (TN)</span>
                    </td>
                </tr>
            </table>
        </div>
        """, unsafe_allow_html=True)

    with col_metrics:
        st.subheader("Performance Metrics")
        
        # Display clinical stats as key-value pairs
        metrics_display = {
            "Accuracy (Concordance)": f"{metrics['Accuracy']*100:.1f}%",
            "Sensitivity (Recall / TPR)": f"{metrics['Sensitivity']*100:.1f}%",
            "Specificity (TNR)": f"{metrics['Specificity']*100:.1f}%",
            "Precision (PPV)": f"{metrics['Precision']*100:.1f}%",
            "Negative Predictive Value (NPV)": f"{metrics['NPV']*100:.1f}%",
            "F1-Score (Harmonic Mean)": f"{metrics['F1_Score']:.3f}",
            "Cohen's Kappa (Agreement)": f"{metrics['Kappa']:.3f}"
        }
        
        for k, v in metrics_display.items():
            st.markdown(f"**{k}:** `{v}`")
            
        st.markdown(f"""
        * **Agreement Strength:** `{get_kappa_interpretation(metrics['Kappa'])}` (Landis & Koch standard).
        """)

    st.markdown("<hr style='border: none; border-top: 1px solid #e1e8ed; margin: 30px 0;' />", unsafe_allow_html=True)

    # Advanced Statistical Testing Section
    st.subheader("Statistical Tests")
    
    # 1. McNemar's Test
    try:
        p_val, method = get_mcnemar_p_value(fp, fn)
        
        # Explain significance
        sig_threshold = 0.05
        significant = p_val < sig_threshold
        
        if significant:
            p_val_desc = f"**p-value:** `{p_val:.5f}` (< 0.05) — 🔴 **Statistically Significant Bias Detected**"
            conclusion = f"The pipeline (PO) shows a significant systematic difference in call rates compared to DO. There is a bias towards either over-calling (FP) or under-calling (FN) deletions."
            alert_func = st.error
        else:
            p_val_desc = f"**p-value:** `{p_val:.5f}` (≥ 0.05) — 🟢 **No Statistically Significant Bias**"
            conclusion = "The errors (False Positives and False Negatives) are relatively balanced. There is no significant evidence of a systematic bias in one direction."
            alert_func = st.success
            
        st.markdown(f"""
        ### McNemar's Test for Paired Nominal Data
        McNemar's test checks whether the pipeline has a systematic bias towards over-calling (False Positives) versus under-calling (False Negatives) deletions.
        - **Null Hypothesis (H0):** The rate of False Positives equals the rate of False Negatives.
        - **Alternative Hypothesis (H1):** The rate of False Positives is different from the rate of False Negatives (systematic bias).
        
        **Test Details:**
        - **Method used:** {method}
        - **Discrepant Calls:** False Positives (FP) = `{fp}`, False Negatives (FN) = `{fn}` (Total = `{fp + fn}`)
        - {p_val_desc}
        """)
        
        alert_func(conclusion)
        
    except Exception as e:
        st.warning(f"Could not calculate McNemar's test. Ensure scipy is correctly configured. Error: {e}")

    # Section: Definitions
    with st.expander("🔬 Metric Definitions & Interpretations for PIs"):
        st.markdown("""
        To assist in validation reports, here are the definitions of the statistical terms used in this portal:
        
        - **Accuracy (Concordance):** The overall percentage of matching results (both Neutral and Deletion) between DO and PO.
        - **Sensitivity (Recall):** The proportion of true deletions (confirmed by DO) that were correctly identified by the pipeline (PO). High sensitivity means very few false negatives (missed deletions).
        - **Specificity:** The proportion of true neutral sites (confirmed by DO) that were correctly identified by the pipeline (PO). High specificity means very few false positives (false alarms).
        - **Precision (PPV):** The proportion of deletions called by the pipeline (PO) that are actually confirmed by DO. High precision means that if the pipeline says "Deletion", it is highly likely to be real.
        - **F1-Score:** The harmonic mean of Precision and Sensitivity. It is a single balanced metric representing the pipeline's effectiveness specifically on deletions.
        - **Cohen's Kappa:** A robust metric that measures the agreement between two clinical tests while correcting for agreement occurring by random chance.
          - `< 0.00`: Poor (Disagreement)
          - `0.00 – 0.20`: Slight agreement
          - `0.21 – 0.40`: Fair agreement
          - `0.41 – 0.60`: Moderate agreement
          - `0.61 – 0.80`: Substantial agreement
          - `0.81 – 1.00`: Almost perfect agreement
        """)

# Run if page file is executed directly (classic fallback)
if __name__ == "__main__":
    run()
elif __name__ == "pages.stats":
    run()
