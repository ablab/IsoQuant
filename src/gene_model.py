
############################################################################
# Copyright (c) 2022-2025 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import json
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import euclidean


def parse_data(data):
    genes = {}
    for condition, condition_data in data.items():
        for gene, gene_data in condition_data.items():
            if gene not in genes:
                genes[gene] = {
                    "chromosome": gene_data["chromosome"],
                    "start": gene_data["start"],
                    "end": gene_data["end"],
                    "strand": gene_data["strand"],
                    "biotype": gene_data["biotype"],
                    "transcripts": {},
                }
            genes[gene]["transcripts"][condition] = gene_data["transcripts"]
            genes[gene][condition] = gene_data["value"]
    return genes


def calculate_deviance(wt_transcripts, condition_transcripts):
    all_transcripts = set(wt_transcripts.keys()).union(
        set(condition_transcripts.keys())
    )

    wt_proportions = [wt_transcripts.get(t, 0) for t in all_transcripts]
    condition_proportions = [condition_transcripts.get(t, 0) for t in all_transcripts]

    total_wt = sum(wt_proportions)
    total_condition = sum(condition_proportions)

    if total_wt > 0:
        wt_proportions = [p / total_wt for p in wt_proportions]
    if total_condition > 0:
        condition_proportions = [p / total_condition for p in condition_proportions]

    distance = euclidean(wt_proportions, condition_proportions)

    # Reduce distance if total unique transcripts are 1
    if len(all_transcripts) == 1:
        distance *= 0.7

    return distance


def calculate_metrics(genes):
    metrics = []
    for gene, gene_data in genes.items():
        wt_transcripts = gene_data["transcripts"].get("wild_type", {})

        for condition in gene_data:
            if condition in [
                "chromosome",
                "start",
                "end",
                "strand",
                "biotype",
                "transcripts",
                "wild_type",
            ]:
                continue
            condition_transcripts = gene_data["transcripts"].get(condition, {})
            deviance = calculate_deviance(wt_transcripts, condition_transcripts)
            metrics.append({"gene": gene, "condition": condition, "deviance": deviance})

            value = gene_data.get(condition, 0)
            wt_value = gene_data.get("wild_type", 0)
            abs_diff = abs(value - wt_value)
            metrics.append(
                {
                    "gene": gene,
                    "condition": condition,
                    "value": value,
                    "abs_diff": abs_diff,
                }
            )

    return pd.DataFrame(metrics)


def check_known_target(gene, known_targets):
    for target in known_targets:
        if "|" in target:
            if any(part in gene for part in target.split("|")):
                return 1
        elif target == gene:
            return 1
    return 0


def rank_genes(df, known_genes_path=None):
    if known_genes_path:
        target_genes_df = pd.read_csv(known_genes_path, header=None, names=["gene"])
        known_targets = target_genes_df["gene"].tolist()
        df["known_target"] = df["gene"].apply(
            lambda x: check_known_target(x, known_targets)
        )
    else:
        df["known_target"] = 0

    value_ranking = df.groupby("gene")["value"].mean().reset_index()
    abs_diff_ranking = df.groupby("gene")["abs_diff"].mean().reset_index()
    deviance_ranking = df.groupby("gene")["deviance"].mean().reset_index()

    value_ranking["rank_value"] = value_ranking["value"].rank(ascending=False)
    abs_diff_ranking["rank_abs_diff"] = abs_diff_ranking["abs_diff"].rank(
        ascending=False
    )
    deviance_ranking["rank_deviance"] = deviance_ranking["deviance"].rank(
        ascending=False
    )

    merged_df = value_ranking[["gene", "rank_value"]].merge(
        abs_diff_ranking[["gene", "rank_abs_diff"]], on="gene"
    )
    merged_df = merged_df.merge(deviance_ranking[["gene", "rank_deviance"]], on="gene")
    merged_df = merged_df.merge(df[["gene", "known_target"]], on="gene")
    # Devalue the importance of overall expression by reducing its weight
    merged_df["combined_rank"] = (
        merged_df["rank_value"]  # Reduced weight for rank_value
        + merged_df["rank_abs_diff"]
        + merged_df["rank_deviance"]
    )

    top_combined_ranking = merged_df.sort_values(by="combined_rank").head(10)
    top_deviance_ranking = merged_df.sort_values(by="rank_deviance").head(10)
    top_100_combined_ranking = merged_df.sort_values(by="combined_rank").head(100)

    return (
        top_combined_ranking,
        top_deviance_ranking,
        top_100_combined_ranking,
        merged_df,
    )


def visualize_ranking(
    top_combined_ranking, top_deviance_ranking, merged_df, output_dir
):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Bar plot for combined rank
    plt.figure(figsize=(12, 8))
    sns.barplot(
        x="combined_rank", y="gene", data=top_combined_ranking, palette="viridis"
    )
    plt.title("Top 10 Genes by Combined Ranking")
    plt.xlabel("Combined Rank")
    plt.ylabel("Gene")
    plt.savefig(os.path.join(output_dir, "top_genes_combined_ranking.png"), dpi=300)
    plt.close()

    # Heatmap for metric ranks
    top_genes = top_combined_ranking["gene"].tolist()
    heatmap_data = merged_df[merged_df["gene"].isin(top_genes)]
    heatmap_data = heatmap_data.set_index("gene")[
        ["rank_value", "rank_abs_diff", "rank_deviance"]
    ]

    plt.figure(figsize=(12, 8))
    sns.heatmap(
        heatmap_data,
        annot=True,
        cmap="RdBu_r",
        linewidths=0.5,
        cbar_kws={"label": "Rank"},
    )
    plt.title("Metric Ranks for Top 10 Genes")
    plt.savefig(os.path.join(output_dir, "metric_ranks_heatmap.png"), dpi=300)
    plt.close()

    # Diverging bar plot for deviance
    plt.figure(figsize=(12, 8))
    sns.barplot(
        x="rank_deviance",
        y="gene",
        data=top_deviance_ranking,
        palette="coolwarm",
        orient="h",
    )
    plt.title("Top 10 Genes by Transcript Deviance from Wild Type")
    plt.xlabel("Rank of Deviance from Wild Type")
    plt.ylabel("Gene")
    plt.axvline(x=0, color="grey", linestyle="--")
    plt.savefig(os.path.join(output_dir, "deviance_from_wild_type.png"), dpi=300)
    plt.close()

    # Scatter plot for rank_value vs rank_abs_diff
    plt.figure(figsize=(12, 8))
    sns.scatterplot(
        x="rank_value",
        y="rank_abs_diff",
        hue="gene",
        data=top_deviance_ranking,
        palette="deep",
        s=100,
    )
    plt.title("Rank Value vs Rank Absolute Difference")
    plt.xlabel("Rank Value")
    plt.ylabel("Rank Absolute Difference")
    plt.savefig(os.path.join(output_dir, "rank_value_vs_rank_abs_diff.png"), dpi=300)
    plt.close()

    # Combined multi-metric visualization
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    sns.barplot(
        x="combined_rank",
        y="gene",
        data=top_combined_ranking,
        palette="viridis",
        ax=axes[0, 0],
    )
    axes[0, 0].set_title("Combined Rank")
    axes[0, 0].set_xlabel("Combined Rank")
    axes[0, 0].set_ylabel("Gene")

    sns.heatmap(
        heatmap_data,
        annot=True,
        cmap="RdBu_r",
        linewidths=0.5,
        cbar_kws={"label": "Rank"},
        ax=axes[0, 1],
    )
    axes[0, 1].set_title("Metric Ranks")

    sns.barplot(
        x="rank_deviance",
        y="gene",
        data=top_deviance_ranking,
        palette="coolwarm",
        orient="h",
        ax=axes[1, 0],
    )
    axes[1, 0].set_title("Transcript Deviance from Wild Type")
    axes[1, 0].set_xlabel("Rank of Deviance from Wild Type")
    axes[1, 0].set_ylabel("Gene")
    axes[1, 0].axvline(x=0, color="grey", linestyle="--")

    sns.scatterplot(
        x="rank_value",
        y="rank_abs_diff",
        hue="gene",
        data=top_deviance_ranking,
        palette="deep",
        s=100,
        ax=axes[1, 1],
    )
    axes[1, 1].set_title("Rank Value vs Rank Absolute Difference")
    axes[1, 1].set_xlabel("Rank Value")
    axes[1, 1].set_ylabel("Rank Absolute Difference")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "combined_visualization.png"), dpi=300)
    plt.close()


def save_top_genes(top_combined_ranking, output_dir, num_genes):
    top_combined_ranking.head(num_genes)[["gene"]].to_csv(
        os.path.join(output_dir, f"top_{num_genes}_genes.txt"),
        index=False,
        header=False,
        sep="\t",
    )
    return os.path.join(output_dir, f"top_{num_genes}_genes.txt")


def rank_and_visualize_genes(
    input_data, output_dir, num_genes=100, known_genes_path=None
):
    genes = parse_data(input_data)
    metrics_df = calculate_metrics(genes)
    top_combined_ranking, top_deviance_ranking, top_100_combined_ranking, merged_df = (
        rank_genes(metrics_df, known_genes_path)
    )
    merged_df = merged_df.drop_duplicates(subset="gene", keep="first")
    top_combined_ranking = merged_df.sort_values(by="combined_rank").head(num_genes)
    top_deviance_ranking = merged_df.sort_values(by="rank_deviance").head(num_genes)

    visualize_ranking(top_combined_ranking, top_deviance_ranking, merged_df, output_dir)
    path = save_top_genes(top_combined_ranking, output_dir, num_genes)

    print(f"\nTop {num_genes} Genes by Combined Ranking:")
    print(top_combined_ranking[["gene", "combined_rank"]])
    print(f"\nDetailed Metrics for Top {num_genes} Genes by Combined Ranking:")
    print(top_combined_ranking)

    merged_df.to_csv(os.path.join(output_dir, "gene_metrics.csv"), index=False)

    return path
