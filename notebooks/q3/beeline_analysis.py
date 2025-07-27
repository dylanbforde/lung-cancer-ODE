import pandas as pd
import numpy as np
import networkx as nx
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


class BeelineAnalyzer:
    def __init__(self, base_path: str = "beeline_outputs/LungCancer_tiny"):
        self.base_path = Path(base_path)
        self.genie3_data = None
        self.ppcor_data = None
        self.genie3_stats = None
        self.ppcor_stats = None

    def load_data(self) -> None:
        # Load GENIE3 data
        genie3_path = self.base_path / "GENIE3" / "outFile.txt"
        if genie3_path.exists():
            self.genie3_data = pd.read_csv(genie3_path, sep="\t")
            print(f"Loaded GENIE3 data: {len(self.genie3_data)} interactions")
        else:
            print(f"GENIE3 file not found: {genie3_path}")

        # Load PPCOR data
        ppcor_path = self.base_path / "PPCOR" / "outFile.txt"
        if ppcor_path.exists():
            self.ppcor_data = pd.read_csv(ppcor_path, sep="\t")
            print(f"Loaded PPCOR data: {len(self.ppcor_data)} interactions")
        else:
            print(f"PPCOR file not found: {ppcor_path}")

    def preprocess_data(self) -> None:
        if self.genie3_data is not None:
            # Sort by importance score
            self.genie3_data = self.genie3_data.sort_values(
                "importance", ascending=False
            )

        if self.ppcor_data is not None:
            # Remove self-correlations and perfect correlations
            self.ppcor_data = self.ppcor_data[
                (self.ppcor_data["Gene1"] != self.ppcor_data["Gene2"])
                & (self.ppcor_data["corVal"] != 1.0)
            ].copy()

            # Sort by absolute correlation value
            self.ppcor_data["abs_corVal"] = abs(self.ppcor_data["corVal"])
            self.ppcor_data = self.ppcor_data.sort_values("abs_corVal", ascending=False)

    def calculate_statistics(self) -> Dict:
        stats = {}

        if self.genie3_data is not None:
            self.genie3_stats = {
                "total_interactions": len(self.genie3_data),
                "unique_genes": len(
                    set(
                        self.genie3_data["TF"].tolist()
                        + self.genie3_data["target"].tolist()
                    )
                ),
                "mean_importance": self.genie3_data["importance"].mean(),
                "std_importance": self.genie3_data["importance"].std(),
                "max_importance": self.genie3_data["importance"].max(),
                "min_importance": self.genie3_data["importance"].min(),
                "top_10_percent_threshold": self.genie3_data["importance"].quantile(
                    0.9
                ),
            }
            stats["GENIE3"] = self.genie3_stats

        if self.ppcor_data is not None:
            self.ppcor_stats = {
                "total_interactions": len(self.ppcor_data),
                "unique_genes": len(
                    set(
                        self.ppcor_data["Gene1"].tolist()
                        + self.ppcor_data["Gene2"].tolist()
                    )
                ),
                "mean_correlation": self.ppcor_data["corVal"].mean(),
                "std_correlation": self.ppcor_data["corVal"].std(),
                "mean_abs_correlation": self.ppcor_data["abs_corVal"].mean(),
                "max_abs_correlation": self.ppcor_data["abs_corVal"].max(),
                "positive_correlations": (self.ppcor_data["corVal"] > 0).sum(),
                "negative_correlations": (self.ppcor_data["corVal"] < 0).sum(),
            }
            stats["PPCOR"] = self.ppcor_stats

        return stats

    def get_top_interactions(self, method: str, n: int = 20) -> pd.DataFrame:
        if method.upper() == "GENIE3" and self.genie3_data is not None:
            return self.genie3_data.head(n)
        elif method.upper() == "PPCOR" and self.ppcor_data is not None:
            return self.ppcor_data.head(n)
        else:
            return pd.DataFrame()

    def find_common_interactions(
        self, genie3_top_n: int = 1000, ppcor_top_n: int = 1000
    ) -> pd.DataFrame:
        if self.genie3_data is None or self.ppcor_data is None:
            return pd.DataFrame()

        # Get top interactions from each method
        genie3_top = self.genie3_data.head(genie3_top_n)
        ppcor_top = self.ppcor_data.head(ppcor_top_n)

        # Create standardized interaction pairs
        genie3_pairs = set()
        for _, row in genie3_top.iterrows():
            pair = tuple(sorted([row["TF"], row["target"]]))
            genie3_pairs.add(pair)

        ppcor_pairs = set()
        for _, row in ppcor_top.iterrows():
            pair = tuple(sorted([row["Gene1"], row["Gene2"]]))
            ppcor_pairs.add(pair)

        # Find common pairs
        common_pairs = genie3_pairs.intersection(ppcor_pairs)

        # Create result dataframe
        common_data = []
        for pair in common_pairs:
            gene1, gene2 = pair

            # Find GENIE3 score
            genie3_score = None
            for _, row in genie3_top.iterrows():
                if (row["TF"] == gene1 and row["target"] == gene2) or (
                    row["TF"] == gene2 and row["target"] == gene1
                ):
                    genie3_score = row["importance"]
                    break

            # Find PPCOR score
            ppcor_score = None
            for _, row in ppcor_top.iterrows():
                if (row["Gene1"] == gene1 and row["Gene2"] == gene2) or (
                    row["Gene1"] == gene2 and row["Gene2"] == gene1
                ):
                    ppcor_score = row["corVal"]
                    break

            common_data.append(
                {
                    "Gene1": gene1,
                    "Gene2": gene2,
                    "GENIE3_importance": genie3_score,
                    "PPCOR_correlation": ppcor_score,
                }
            )

        return pd.DataFrame(common_data)

    def get_gene_degree_distribution(self, method: str, top_n: int = 1000) -> pd.Series:
        if method.upper() == "GENIE3" and self.genie3_data is not None:
            top_data = self.genie3_data.head(top_n)
            genes = list(top_data["TF"]) + list(top_data["target"])
        elif method.upper() == "PPCOR" and self.ppcor_data is not None:
            top_data = self.ppcor_data.head(top_n)
            genes = list(top_data["Gene1"]) + list(top_data["Gene2"])
        else:
            return pd.Series()

        degree_counts = pd.Series(genes).value_counts()
        return degree_counts

    def create_network_graph(self, method: str, top_n: int = 100) -> nx.Graph:
        G = nx.Graph()

        if method.upper() == "GENIE3" and self.genie3_data is not None:
            top_data = self.genie3_data.head(top_n)
            for _, row in top_data.iterrows():
                G.add_edge(row["TF"], row["target"], weight=row["importance"])
        elif method.upper() == "PPCOR" and self.ppcor_data is not None:
            top_data = self.ppcor_data.head(top_n)
            for _, row in top_data.iterrows():
                G.add_edge(row["Gene1"], row["Gene2"], weight=abs(row["corVal"]))

        return G

    def plot_score_distributions(
        self, figsize: Tuple[int, int] = (12, 5)
    ) -> plt.Figure:
        fig, axes = plt.subplots(1, 2, figsize=figsize)

        if self.genie3_data is not None:
            axes[0].hist(
                self.genie3_data["importance"], bins=50, alpha=0.7, color="blue"
            )
            axes[0].set_title("GENIE3 Importance Score Distribution")
            axes[0].set_xlabel("Importance Score")
            axes[0].set_ylabel("Frequency")
            axes[0].axvline(
                self.genie3_stats["mean_importance"],
                color="red",
                linestyle="--",
                label=f"Mean: {self.genie3_stats['mean_importance']:.6f}",
            )
            axes[0].legend()

        if self.ppcor_data is not None:
            axes[1].hist(self.ppcor_data["corVal"], bins=50, alpha=0.7, color="green")
            axes[1].set_title("PPCOR Correlation Distribution")
            axes[1].set_xlabel("Correlation Value")
            axes[1].set_ylabel("Frequency")
            axes[1].axvline(
                self.ppcor_stats["mean_correlation"],
                color="red",
                linestyle="--",
                label=f"Mean: {self.ppcor_stats['mean_correlation']:.6f}",
            )
            axes[1].legend()

        plt.tight_layout()
        return fig

    def plot_top_genes_degree(
        self, method: str, top_n: int = 20, figsize: Tuple[int, int] = (10, 6)
    ) -> plt.Figure:
        """Plot degree distribution for top genes"""
        degree_dist = self.get_gene_degree_distribution(method, top_n=1000)
        top_genes = degree_dist.head(top_n)

        fig, ax = plt.subplots(figsize=figsize)
        top_genes.plot(kind="bar", ax=ax)
        ax.set_title(f"Top {top_n} Genes by Degree ({method})")
        ax.set_xlabel("Genes")
        ax.set_ylabel("Degree (Number of Connections)")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        return fig

    def get_top_genes_for_survival(
        self, method: str = "GENIE3", top_n: int = 50
    ) -> List[str]:
        """Extract top genes for survival modeling based on network centrality"""
        if method.upper() == "GENIE3" and self.genie3_data is not None:
            degree_dist = self.get_gene_degree_distribution("GENIE3", top_n=1000)
            top_genes = degree_dist.head(top_n).index.tolist()
        elif method.upper() == "PPCOR" and self.ppcor_data is not None:
            degree_dist = self.get_gene_degree_distribution("PPCOR", top_n=1000)
            top_genes = degree_dist.head(top_n).index.tolist()
        else:
            return []

        return top_genes

    def get_high_importance_genes(
        self, method: str = "GENIE3", threshold_percentile: float = 0.95
    ) -> List[str]:
        """Get genes involved in high-importance interactions"""
        genes = set()

        if method.upper() == "GENIE3" and self.genie3_data is not None:
            threshold = self.genie3_data["importance"].quantile(threshold_percentile)
            high_imp = self.genie3_data[self.genie3_data["importance"] >= threshold]
            genes.update(high_imp["TF"].tolist())
            genes.update(high_imp["target"].tolist())
        elif method.upper() == "PPCOR" and self.ppcor_data is not None:
            threshold = self.ppcor_data["abs_corVal"].quantile(threshold_percentile)
            high_cor = self.ppcor_data[self.ppcor_data["abs_corVal"] >= threshold]
            genes.update(high_cor["Gene1"].tolist())
            genes.update(high_cor["Gene2"].tolist())

        return list(genes)

    def save_visualizations(self, output_dir: str = "visualizations") -> None:
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

        # Score distributions
        fig = self.plot_score_distributions()
        fig.savefig(
            output_path / "score_distributions.png", dpi=300, bbox_inches="tight"
        )
        plt.close(fig)

        # Top genes degree plots
        if self.genie3_data is not None:
            fig = self.plot_top_genes_degree("GENIE3")
            fig.savefig(
                output_path / "genie3_top_genes_degree.png",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close(fig)

        if self.ppcor_data is not None:
            fig = self.plot_top_genes_degree("PPCOR")
            fig.savefig(
                output_path / "ppcor_top_genes_degree.png", dpi=300, bbox_inches="tight"
            )
            plt.close(fig)

        # Network comparison plot
        fig = self.plot_method_comparison()
        fig.savefig(output_path / "method_comparison.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

        # Common interactions heatmap
        if self.genie3_data is not None and self.ppcor_data is not None:
            fig = self.plot_common_interactions_heatmap()
            fig.savefig(
                output_path / "common_interactions_heatmap.png",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close(fig)

        print(f"Visualizations saved to {output_path}")

    def plot_method_comparison(self, figsize: Tuple[int, int] = (12, 8)) -> plt.Figure:
        fig, axes = plt.subplots(2, 2, figsize=figsize)

        if self.genie3_data is not None:
            # GENIE3 score vs rank
            genie3_scores = self.genie3_data["importance"].values
            axes[0, 0].plot(range(len(genie3_scores)), genie3_scores, "b-", alpha=0.7)
            axes[0, 0].set_title("GENIE3: Importance Score vs Rank")
            axes[0, 0].set_xlabel("Rank")
            axes[0, 0].set_ylabel("Importance Score")
            axes[0, 0].set_yscale("log")

            # GENIE3 degree distribution
            degree_dist = self.get_gene_degree_distribution("GENIE3", 1000)
            axes[1, 0].hist(degree_dist.values, bins=30, alpha=0.7, color="blue")
            axes[1, 0].set_title("GENIE3: Gene Degree Distribution")
            axes[1, 0].set_xlabel("Degree")
            axes[1, 0].set_ylabel("Number of Genes")

        if self.ppcor_data is not None:
            # PPCOR score vs rank
            ppcor_scores = self.ppcor_data["abs_corVal"].values
            axes[0, 1].plot(range(len(ppcor_scores)), ppcor_scores, "g-", alpha=0.7)
            axes[0, 1].set_title("PPCOR: |Correlation| vs Rank")
            axes[0, 1].set_xlabel("Rank")
            axes[0, 1].set_ylabel("Absolute Correlation")
            axes[0, 1].set_yscale("log")

            # PPCOR degree distribution
            degree_dist = self.get_gene_degree_distribution("PPCOR", 1000)
            axes[1, 1].hist(degree_dist.values, bins=30, alpha=0.7, color="green")
            axes[1, 1].set_title("PPCOR: Gene Degree Distribution")
            axes[1, 1].set_xlabel("Degree")
            axes[1, 1].set_ylabel("Number of Genes")

        plt.tight_layout()
        return fig

    def plot_common_interactions_heatmap(
        self, top_n: int = 500, figsize: Tuple[int, int] = (10, 8)
    ) -> plt.Figure:
        common_interactions = self.find_common_interactions(top_n, top_n)

        if len(common_interactions) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(
                0.5,
                0.5,
                "No common interactions found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title(f"Common Interactions Heatmap (Top {top_n})")
            return fig

        # Create correlation matrix for visualization
        genes = list(
            set(
                common_interactions["Gene1"].tolist()
                + common_interactions["Gene2"].tolist()
            )
        )

        if len(genes) > 50:  # Limit to prevent overcrowded plot
            genes = genes[:50]

        # Create adjacency matrix
        adj_matrix = pd.DataFrame(0, index=genes, columns=genes)

        for _, row in common_interactions.iterrows():
            if row["Gene1"] in genes and row["Gene2"] in genes:
                adj_matrix.loc[row["Gene1"], row["Gene2"]] = row["GENIE3_importance"]
                adj_matrix.loc[row["Gene2"], row["Gene1"]] = row["GENIE3_importance"]

        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(adj_matrix, annot=False, cmap="viridis", ax=ax)
        ax.set_title(f"Common Interactions Network (Top {top_n})")
        plt.tight_layout()
        return fig

    def generate_summary_report(self) -> str:
        report = "=== BEELINE ANALYSIS SUMMARY ===\n\n"

        if self.genie3_stats:
            report += "GENIE3 Results:\n"
            report += (
                f"  - Total interactions: {self.genie3_stats['total_interactions']:,}\n"
            )
            report += f"  - Unique genes: {self.genie3_stats['unique_genes']}\n"
            report += (
                f"  - Mean importance: {self.genie3_stats['mean_importance']:.6f}\n"
            )
            report += f"  - Max importance: {self.genie3_stats['max_importance']:.6f}\n"
            report += f"  - Top 10% threshold: {self.genie3_stats['top_10_percent_threshold']:.6f}\n\n"

        if self.ppcor_stats:
            report += "PPCOR Results:\n"
            report += (
                f"  - Total interactions: {self.ppcor_stats['total_interactions']:,}\n"
            )
            report += f"  - Unique genes: {self.ppcor_stats['unique_genes']}\n"
            report += (
                f"  - Mean correlation: {self.ppcor_stats['mean_correlation']:.6f}\n"
            )
            report += f"  - Mean absolute correlation: {self.ppcor_stats['mean_abs_correlation']:.6f}\n"
            report += f"  - Positive correlations: {self.ppcor_stats['positive_correlations']:,}\n"
            report += f"  - Negative correlations: {self.ppcor_stats['negative_correlations']:,}\n\n"

        # Add top interactions
        if self.genie3_data is not None:
            report += "Top 10 GENIE3 Interactions:\n"
            top_genie3 = self.get_top_interactions("GENIE3", 10)
            for _, row in top_genie3.iterrows():
                report += f"  {row['TF']} -> {row['target']}: {row['importance']:.6f}\n"
            report += "\n"

        if self.ppcor_data is not None:
            report += "Top 10 PPCOR Interactions:\n"
            top_ppcor = self.get_top_interactions("PPCOR", 10)
            for _, row in top_ppcor.iterrows():
                report += f"  {row['Gene1']} <-> {row['Gene2']}: {row['corVal']:.6f}\n"
            report += "\n"

        return report

    def save_gene_lists_for_survival(
        self, output_dir: str = "gene_lists"
    ) -> Dict[str, List[str]]:
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

        gene_lists = {}

        # Top hub genes (high degree)
        if self.genie3_data is not None:
            genie3_hubs = self.get_top_genes_for_survival("GENIE3", 50)
            gene_lists["genie3_hub_genes"] = genie3_hubs
            pd.Series(genie3_hubs).to_csv(
                output_path / "genie3_hub_genes.csv", header=["gene"], index=False
            )

        if self.ppcor_data is not None:
            ppcor_hubs = self.get_top_genes_for_survival("PPCOR", 50)
            gene_lists["ppcor_hub_genes"] = ppcor_hubs
            pd.Series(ppcor_hubs).to_csv(
                output_path / "ppcor_hub_genes.csv", header=["gene"], index=False
            )

        # High importance/correlation genes
        if self.genie3_data is not None:
            genie3_high = self.get_high_importance_genes("GENIE3", 0.95)
            gene_lists["genie3_high_importance"] = genie3_high
            pd.Series(genie3_high).to_csv(
                output_path / "genie3_high_importance_genes.csv",
                header=["gene"],
                index=False,
            )

        if self.ppcor_data is not None:
            ppcor_high = self.get_high_importance_genes("PPCOR", 0.95)
            gene_lists["ppcor_high_correlation"] = ppcor_high
            pd.Series(ppcor_high).to_csv(
                output_path / "ppcor_high_correlation_genes.csv",
                header=["gene"],
                index=False,
            )

        # Common genes between methods
        if self.genie3_data is not None and self.ppcor_data is not None:
            common_hubs = list(
                set(gene_lists.get("genie3_hub_genes", []))
                & set(gene_lists.get("ppcor_hub_genes", []))
            )
            gene_lists["common_hub_genes"] = common_hubs
            pd.Series(common_hubs).to_csv(
                output_path / "common_hub_genes.csv", header=["gene"], index=False
            )

        print(f"Gene lists saved to {output_path}")
        return gene_lists


def main():
    analyzer = BeelineAnalyzer()
    analyzer.load_data()
    analyzer.preprocess_data()
    stats = analyzer.calculate_statistics()

    print(analyzer.generate_summary_report())

    # Find common interactions
    common = analyzer.find_common_interactions(1000, 1000)
    print(f"Common interactions in top 1000: {len(common)}")

    return analyzer


if __name__ == "__main__":
    analyzer = main()
