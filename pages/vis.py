import streamlit as st
from utils.plotting import plot_cohort_heatmap, plot_class_distribution, plot_sensitivity_specificity
from utils.statistics import get_gene_metrics

def run():
    st.markdown("""
    <div class="portal-header">
        <h1>📈 Interactive Visualization Suite</h1>
        <p>Explore cohort patterns, deletion co-occurrences, and validation metrics interactively.</p>
    </div>
    """, unsafe_allow_html=True)

    long_df = st.session_state["long_df"]
    genes = st.session_state["genes"]
    
    # Gene metrics needed for sensitivity/specificity plot
    gene_metrics_df = get_gene_metrics(long_df)

    # Tabs for different views
    tab_heatmap, tab_breakdown, tab_tradeoffs = st.tabs([
        "🔬 Clinical Cohort Heatmap", 
        "📊 Comparison Breakdown by Gene", 
        "🎯 Sensitivity vs. Specificity Space"
    ])

    with tab_heatmap:
        st.markdown("""
        ### Clinical Cohort Heatmap
        This heatmap represents the classification output across the entire cohort (patients as rows, genes as columns).
        Use the controls below to sort patients or filter the cohort view.
        """)
        
        col_sort, col_filter = st.columns([1, 2])
        with col_sort:
            sort_by = st.selectbox(
                "Sort Cohort (Rows) By",
                ["RegId", "Concordance", "Number of Deletions"],
                index=0,
                key="heatmap_sort"
            )
        with col_filter:
            selected_heatmap_genes = st.multiselect(
                "Filter Genes displayed in Heatmap",
                genes,
                default=genes,
                key="heatmap_genes"
            )
            
        if selected_heatmap_genes:
            fig_heatmap = plot_cohort_heatmap(long_df, selected_heatmap_genes, sort_by=sort_by)
            st.plotly_chart(fig_heatmap, width="stretch")
        else:
            st.warning("Please select at least one gene to display the heatmap.")

    with tab_breakdown:
        st.markdown("""
        ### Comparison Breakdown per Gene
        This stacked bar chart shows the proportion of matching and discordant calls for each gene.
        """)
        
        selected_breakdown_genes = st.multiselect(
            "Filter Genes to display",
            genes,
            default=genes,
            key="breakdown_genes"
        )
        
        if selected_breakdown_genes:
            fig_breakdown = plot_class_distribution(long_df, selected_breakdown_genes)
            st.plotly_chart(fig_breakdown, width="stretch")
        else:
            st.warning("Please select at least one gene to display.")

    with tab_tradeoffs:
        st.markdown("""
        ### Gene Validation Space: Sensitivity vs. Specificity
        This multi-dimensional scatter plot displays each gene as a data point.
        - **X-axis:** Sensitivity (True Deletion Detection Rate)
        - **Y-axis:** Specificity (True Neutral Detection Rate)
        - **Marker Color:** Cohen's Kappa score (agreement corrected for chance)
        - **Marker Size:** F1-score (overall metric for deletions)
        
        *Perfect agreement is at the top-right corner (1.0, 1.0).*
        """)
        
        fig_tradeoffs = plot_sensitivity_specificity(gene_metrics_df)
        st.plotly_chart(fig_tradeoffs, width="stretch")

# Run if page file is executed directly (classic fallback)
if __name__ == "__main__":
    run()
elif __name__ == "pages.vis":
    run()
