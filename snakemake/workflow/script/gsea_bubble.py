"""
example: python gsea_bubble.py -edbs  \
gsea/Comparison1.rnk.txt/m5.go.bp.v2024.1.Mm.symbols.gmt.GseaPreranked/edb/results.edb \
gsea/Comparison2.rnk.txt/m5.go.bp.v2024.1.Mm.symbols.gmt.GseaPreranked/edb/results.edb \
-output gsea_bubble/go.bp.pdf \
-alpha 0.05 \
-topn 100
"""

import numpy as np
import pandas as pd
import seaborn as sns
import re
import os
import math
import matplotlib.pyplot as plt
import argparse


def parse_gsea_edb(edb_path):
    """
    Reads a GSEA .edb file and extracts relevant information.
    Input example:
    <DTG RANKED_LIST="Comparison1.rnk" TEMPLATE="na_as_pre_ranked" GENESET="gene_sets.gmt#CHR19A" ES="-0.1290" \
    NES="-2.2559" NP="0.0000" FDR="0.0058" FWER="0.0030" RND_ES="0.0601 -0.0291 0.0570 -0.0666

    Args:
        edb_path (str): Path to the GSEA .edb file.

    Returns:
        pd.DataFrame: DataFrame containing the extracted values.
    """
    data = []

    with open(edb_path, "r") as file:
        for line in file:
            match = re.search(
                r'RANKED_LIST="(.+?)"\s+TEMPLATE=".+?"\s+GENESET=".+?#(.+?)"\s+ES="([-0-9.]+)"\s+NES="([-0-9.]+)"\s+NP="([-0-9.]+)"\s+FDR="([-0-9.]+)"',
                line
            )
            if match:
                ranked_list, geneset, es, nes, np, fdr = match.groups()
                # if not (math.isnan(es) or math.isnan(nes) or math.isnan(np) or math.isnan(fdr)):
                if not (math.isnan(float(nes)) or math.isnan(float(fdr))):  # skip row if NES or FDR is NaN
                    ranked_list = ranked_list.replace('.rnk', '')
                    fdr = float(fdr) if fdr else 1.0  # Default FDR to 1.0 if not present
                    np = float(np) if np else 1.0
                    data.append([ranked_list, geneset, float(es), float(nes), np, fdr])
                else:
                    print(f"Skipping row with NaN values: {line.strip()}")

    df = pd.DataFrame(data, columns=["Comparison", "GeneSet", "ES", "NES", "NP", "FDR"])
    return df


def filter_top_gsea_results(df, alpha=0.05, topn=100):
    """
    Filters GSEA results to keep only GeneSets where at least one comparison has FDR < alpha.
    If more than `topn` GeneSets pass, select the top `topn` by minimum FDR.

    Args:
        df (pd.DataFrame): Melted dataframe with columns ['Comparison', 'GeneSet', 'NES', 'FDR', 'p-value'].
        alpha (float): Threshold for filtering GeneSets with FDR (default 0.05).
        topn (int): Maximum number of GeneSets to keep (default 100).

    Returns:
        pd.DataFrame: Filtered dataframe.
    """
    n_comparisons = len(df['Comparison'].unique())

    # Step 1: Find GeneSets where at least one comparison has FDR < alpha
    significant_gene_sets = df[df["FDR"] < alpha]["GeneSet"].unique()
    df_filtered = df[df["GeneSet"].isin(significant_gene_sets)]

    # Step 2: Filter out GeneSets with result rows less than n_comparisons
    df_filtered = df_filtered.groupby("GeneSet").filter(lambda x: len(x) >= n_comparisons)

    # Step 3: If more than topn GeneSets, rank by min FDR and keep topn
    if len(significant_gene_sets) > topn:
        # Compute min FDR for each GeneSet across all comparisons
        gene_set_min_fdr = df_filtered.groupby("GeneSet")["FDR"].min().reset_index()

        # Select the topn GeneSets with the smallest min FDR
        top_gene_sets = gene_set_min_fdr.nsmallest(topn, "FDR")["GeneSet"]

        # Filter the original dataframe to keep only these top GeneSets
        df_filtered = df_filtered[df_filtered["GeneSet"].isin(top_gene_sets)]

    return df_filtered.sort_values(["GeneSet", "Comparison"])


def create_bubble_plot(df, output_path="folder/plot.pdf", alpha='alpha'):
    """
    Create a multi-bubble plot from GSEA results.
    Args:
        df (pd.DataFrame): DataFrame containing GSEA results with \
        columns ['comparison_name', 'GSEA_DB', 'NES', 'p_value'].
        output_path (str): Path to save the plot.
        alpha (float): Significance level for filtering GeneSets. Not used for filtering, just for the title.

    Returns:
        None
    """
    # Create bubble plot
    nrows = df['GeneSet'].unique().shape[0]
    plt.figure(figsize=(8, nrows * 0.2 + 2))

    sns.scatterplot(
        data=df,
        x="NES",
        y="GeneSet",
        size="-log10(FDR)",
        hue="Comparison",
        palette="Set1",
        sizes=(20, 80),  # Control bubble size range
        edgecolor="black",
        alpha=0.7
    )

    # Add vertical reference lines
    plt.axvline(x=0, linestyle="--", color="grey", alpha=0.6)  # NES=0 line
    plt.axvline(x=-3, linestyle="--", color="grey", alpha=0.6)
    plt.axvline(x=3, linestyle="--", color="grey", alpha=0.6)  # NES=3 line, confident range

    # Improve grid visibility
    plt.grid(axis='y', linestyle='-', linewidth=0.6, alpha=0.7)  # Light horizontal grid
    plt.gca().set_axisbelow(True)  # Ensure grid is behind the scatter points

    # Set x-axis limits
    MAX_X_RANGE = df['NES'].abs().max() + 0.5
    print(MAX_X_RANGE)
    plt.xlim(-MAX_X_RANGE, MAX_X_RANGE)

    plt.xlabel("Normalized Enrichment Score (NES)")
    plt.ylabel("Gene Set")
    plt.title("Multi-Bubble Plot of GSEA Results\n" +
              f"Top {len(df['GeneSet'].unique())} Gene Sets with min-FDR < {alpha}")

    # Remove default size legend as they are not integers
    filtered_handles = []
    filtered_labels = []
    handles, labels = plt.gca().get_legend_handles_labels()
    for h, l in zip(handles, labels):
        if l == "-log10(FDR)":
            break
        filtered_handles.append(h)
        filtered_labels.append(l)

    # Manually create size legend for -log10(FDR)
    size_handles = [
        plt.Line2D([0], [0], color="white", linestyle=""),  # Dummy for spacing
        plt.Line2D([0], [0], color="white", linestyle=""),  # Dummy for spacing
        plt.scatter([], [], s=20, color='gray', edgecolor='black', label='1'),
        plt.scatter([], [], s=40, color='gray', edgecolor='black', label='2'),
        plt.scatter([], [], s=60, color='gray', edgecolor='black', label='3'),
        plt.scatter([], [], s=80, color='gray', edgecolor='black', label='4')
    ]
    size_labels = [" ", "-log10(FDR)", "1", "2", "3", "4"]

    merged_handles = filtered_handles + size_handles
    merged_labels = filtered_labels + size_labels

    plt.legend(merged_handles, merged_labels, bbox_to_anchor=(1, 1))

    plt.savefig(output_path, format="pdf", bbox_inches="tight")


def main():
    parser = argparse.ArgumentParser(description="Parse GSEA .edb files and create a bubble plot.")
    parser.add_argument("-edbs", nargs="+", required=True, help="List of .edb files.")
    parser.add_argument("-output", required=True, help="Output fname for the plot.")
    parser.add_argument("-alpha", type=float, default=0.05, help="max-value for: min-FDR in all comparisons; threshold for filtering GeneSets (default: 0.05).")
    parser.add_argument("-topn", type=int, default=1000, help="Top n GeneSets to keep based on ranked min-FDR (default: 1000).")
    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Read all .edb files
    all_data = [parse_gsea_edb(edb) for edb in args.edbs]
    combined_df = pd.concat(all_data, ignore_index=True)
    combined_df["-log10(FDR)"] = -np.log10(combined_df["FDR"].replace(0, 1e-4))  # Replace 0 with a small value to avoid log(0) issue
    df_filtered = filter_top_gsea_results(combined_df, alpha=args.alpha, topn=args.topn)

    # Create the plot
    create_bubble_plot(df_filtered, args.output, args.alpha)

if __name__ == "__main__":
    main()
