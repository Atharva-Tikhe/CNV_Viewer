import pandas as pd
import numpy as np

def calculate_metrics(tp, tn, fp, fn):
    """
    Calculates standard performance metrics based on TP, TN, FP, FN.
    """
    total = tp + tn + fp + fn
    accuracy = (tp + tn) / total if total > 0 else 0
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    
    f1 = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
    
    # Cohen's Kappa
    # Observed agreement (Po)
    po = (tp + tn) / total if total > 0 else 0
    # Expected agreement (Pe)
    if total > 0:
        p_yes = ((tp + fp) / total) * ((tp + fn) / total)
        p_no = ((fn + tn) / total) * ((fp + tn) / total)
        pe = p_yes + p_no
        kappa = (po - pe) / (1 - pe) if pe < 1 else 1.0
    else:
        kappa = 0
        
    return {
        "TP": tp,
        "TN": tn,
        "FP": fp,
        "FN": fn,
        "Total": total,
        "Accuracy": accuracy,
        "Sensitivity": sensitivity,
        "Specificity": specificity,
        "Precision": precision,
        "NPV": npv,
        "F1_Score": f1,
        "Kappa": kappa
    }

def get_kappa_interpretation(kappa):
    """
    Returns Landis & Koch interpretation of Cohen's Kappa.
    """
    if kappa < 0:
        return "Poor (Disagreement)"
    elif kappa <= 0.20:
        return "Slight"
    elif kappa <= 0.40:
        return "Fair"
    elif kappa <= 0.60:
        return "Moderate"
    elif kappa <= 0.80:
        return "Substantial"
    else:
        return "Almost Perfect"

def get_gene_metrics(long_df):
    """
    Calculates metrics for each gene separately.
    Returns a pandas DataFrame sorted by Gene.
    """
    genes = sorted(long_df["Gene"].unique())
    gene_metrics_list = []
    
    for gene in genes:
        gene_data = long_df[long_df["Gene"] == gene]
        tp = (gene_data["Class"] == "TP").sum()
        tn = (gene_data["Class"] == "TN").sum()
        fp = (gene_data["Class"] == "FP").sum()
        fn = (gene_data["Class"] == "FN").sum()
        
        metrics = calculate_metrics(tp, tn, fp, fn)
        metrics["Gene"] = gene
        gene_metrics_list.append(metrics)
        
    res_df = pd.DataFrame(gene_metrics_list)
    # Reorder columns
    cols = ["Gene", "TP", "TN", "FP", "FN", "Total", "Accuracy", "Sensitivity", "Specificity", "Precision", "NPV", "F1_Score", "Kappa"]
    return res_df[cols]

def get_overall_metrics(long_df):
    """
    Calculates overall metrics aggregated across all genes.
    """
    tp = (long_df["Class"] == "TP").sum()
    tn = (long_df["Class"] == "TN").sum()
    fp = (long_df["Class"] == "FP").sum()
    fn = (long_df["Class"] == "FN").sum()
    
    return calculate_metrics(tp, tn, fp, fn)
