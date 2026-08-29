import pandas as pd
import os

panel = [
    "BTG1",
    "CDKN2A",
    "CDKN2B",
    "CDKN2A/B",
    "EBF1",
    "ETV6",
    "IKZF1",
    "PAX5",
    "RB1",
    "PAR1",
]

pipeline_output = [
    "BTG1",
    "CDKN2A",
    "CDKN2B",
    "EBF1",
    "ETV6",
    "IKZF1",
    "PAX5",
    "RB1",
    "PAR1",
]

PAR_genes = ["SHOX", "CRLF2", "IL3RA", "ASMTL", "P2RY8"]

# update panel, pipeline genes from original panel
df = pd.read_csv("/home/atharva/opt/panel_genes.bed", sep="\t")

panel_genes = list(df["gene"])

for gene in panel_genes:
    if gene not in panel and gene not in PAR_genes:
        panel.append(gene)
        pipeline_output.append(gene)


def get_cn_results(regid, df):
    ukall = [regid]
    genes = {}

    for idx, row in df[df["RegId"] == regid].iterrows():
        if row["Target"] == "Copy Number Profile":
            ukall.append((row["Gene"], row["Result"]))

        if row["Gene"] in panel and row["Target"] == "Deletion":
            genes[row["Gene"]] = "D"

        elif row["Gene"] in panel and row["Target"] == "Amplification":
            genes[row["Gene"]] = "A"

    for gene in panel:
        if gene not in list(genes.keys()):
            genes[gene] = "N"

    if genes["CDKN2A/B"] == "D":
        genes["CDKN2A"] = "D"
        genes["CDKN2B"] = "D"

    return ukall, genes


def get_pipeline_results(sample, executions_path):
    try:
        pipeline_results = pd.read_csv(
            f"{executions_path}/{sample}/stats/gene_status.csv", header=0
        )
        pipeline_res = {}

        for gene in pipeline_output:
            if gene in list(pipeline_results.columns):
                pipeline_res[gene] = (
                    "D" if pipeline_results.loc[:, gene].values[0] == "Yes" else "N"
                )
            else:
                pipeline_res[gene] = "N"

        # The pipeline treats CDKN2A/B separately but the reports mention both together, so we merge same results
        if pipeline_res["CDKN2A"] == pipeline_res["CDKN2B"]:
            pipeline_res["CDKN2A/B"] = pipeline_res["CDKN2A"]
        else:
            pipeline_res["CDKN2A/B"] = "N"

        return pipeline_res
    except Exception:
        print(f"{sample} is incomplete")


def matches(diag_out, exec_path, out_path):
    do = diag_out
    executions_path = exec_path
    output_path = out_path

    df = pd.read_csv(do)
    df["RegId"] = df["RegId"].astype(int)

    executions = os.listdir(executions_path)
    if "cnv_data.csv" in executions:
        executions.remove("cnv_data.csv")

    records = []
    for sample in executions:
        sid = sample.split("_")[1]
        sid = int(sid)

        ukall, genes = get_cn_results(sid, df)
        pipeline_res = get_pipeline_results(sample, executions_path)

        gdf = pd.DataFrame(genes, index=["DO"])
        gdf["RegId"] = sid
        pdf = pd.DataFrame(pipeline_res, index=["PO"])
        pdf["RegId"] = sid

        records.append(gdf)
        records.append(pdf)

        print(f"processed {sample}")

    dfs = pd.concat(records)
    dfs.to_csv(output_path)


# matches()
