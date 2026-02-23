import logging
from typing import List, TypedDict, Union

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.sparse import issparse

# --- Logger Configuration ---
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# --- Global Constants ---
# Gene lists for cell cycle scoring, based on Tirosh et al., 2015.
S_GENES = [
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM2",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDCA7",
    "DTL",
    "PRIM1",
    "UHRF1",
    "MLF1IP",
    "HELLS",
    "RFC2",
    "RPA2",
    "NASP",
    "RAD51AP1",
    "GMNN",
    "WDR76",
    "SLBP",
    "CCNE2",
    "UBR7",
    "POLD3",
    "MSH2",
    "ATAD2",
    "RAD51",
    "RRM2",
    "CDC45",
    "CDC6",
    "EXO1",
    "TIPIN",
    "DSCC1",
    "BLM",
    "CASP8AP2",
    "USP1",
    "CLSPN",
    "POLA1",
    "CHAF1B",
    "BRIP1",
    "E2F8",
]

G2M_GENES = [
    "HMGB2",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "NUF2",
    "CKS1B",
    "MKI67",
    "TMPO",
    "CENPF",
    "TACC3",
    "FAM64A",
    "SMC4",
    "CCNB2",
    "CKAP2L",
    "CKAP2",
    "AURKB",
    "BUB1",
    "KIF11",
    "ANP32E",
    "TUBB4B",
    "GTSE1",
    "KIF20B",
    "HJURP",
    "CDCA3",
    "HN1",
    "CDC20",
    "TTK",
    "CDC25C",
    "KIF2C",
    "RANGAP1",
    "NCAPD2",
    "DLGAP5",
    "CDCA2",
    "CDCA8",
    "ECT2",
    "KIF23",
    "HMMR",
    "AURKA",
    "PSRC1",
    "ANLN",
    "LBR",
    "CKAP5",
    "CENPE",
    "CTCF",
    "NEK2",
    "G2E3",
    "GAS2L3",
    "CBX5",
    "CENPA",
]


# --- Custom Type Definitions ---
# Using TypedDict to define the structure of the return value for analyze_expression.
class ExpressionResults(TypedDict):
    """A dictionary containing the results of the expression analysis."""

    individual_values: pd.DataFrame
    averaged_values: pd.DataFrame
    percent_above_median: pd.DataFrame


# --- Helper Functions ---
def _ensure_in(container, items: List[str], name: str):
    """
    A private helper function to check if all items are in a given container.
    Raises a ValueError with a formatted message if any items are missing.
    This promotes code reuse and keeps validation logic separate.

    Args:
        container: The collection to check against (e.g., list, DataFrame columns).
        items: A list of strings that should be in the container.
        name: A descriptive name for the items being checked (for the error message).
    """
    missing = [i for i in items if i not in container]
    if missing:
        raise ValueError(f"{name} not found: {', '.join(missing)}")


# --- Main Analysis Functions ---
def analyze_expression(
    adata: ad.AnnData,
    obs_keys: Union[str, List[str]],
    var_names: Union[str, List[str]],
    plot: bool = True,
    plot_type: str = "violin",
    save_dir: Union[str, None] = None
) -> ExpressionResults:
    """
    Analyzes and optionally visualizes gene expression from an AnnData object.

    Args:
        adata: The annotated data matrix.
        obs_keys: Key(s) from adata.obs to group cells by.
        var_names: Gene name(s) from adata.var_names to analyze.
        plot: If True, generates a plot for each gene.
        plot_type: The type of plot to generate ('violin' or 'box').

    Returns:
        An ExpressionResults dictionary containing detailed and summary DataFrames.
    """
    # --- 1. Input Validation and Standardization ---
    # Ensure obs_keys and var_names are lists for consistent processing.
    if isinstance(obs_keys, str):
        obs_keys = [obs_keys]
    if isinstance(var_names, str):
        var_names = [var_names]

    # Use the helper function to perform validation.
    _ensure_in(adata.obs.columns, obs_keys, "Observation key(s)")
    _ensure_in(adata.var_names, var_names, "Variable(s)")

    logger.info(
        f"Analyzing {len(var_names)} variable(s) grouped by {', '.join(obs_keys)}..."
    )

    # --- 2. Data Extraction ---
    # Create the base DataFrame with observation metadata.
    df = adata.obs[obs_keys].copy()
    # Subset the AnnData object to only the genes of interest for efficiency.
    sub_adata = adata[:, var_names]

    # Helper function to safely get data from a layer or default to .X.
    def _get_layer_or_X(layer_name: str):
        # Default to .X if the specified layer doesn't exist.
        layer = sub_adata.layers.get(layer_name, sub_adata.X)
        return layer.toarray() if issparse(layer) else layer

    # Extract log-normalized data, assuming it's in a layer 'l1p' or in .X.
    lognorm = _get_layer_or_X("l1p")
    for i, var in enumerate(sub_adata.var_names):
        df[f"{var}_lognorm"] = lognorm[:, i]

    # Safely extract raw counts from the 'counts' layer if it exists.
    raw = sub_adata.layers.get("counts")
    if raw is not None:
        raw = raw.toarray() if issparse(raw) else raw
        for i, var in enumerate(sub_adata.var_names):
            df[f"{var}_raw"] = raw[:, i]
    else:
        logger.warning("Raw counts layer 'counts' not found.")
        return "Raw counts layer 'counts' not found"

    # --- 3. Perform Calculations ---
    lognorm_cols = [f"{v}_lognorm" for v in var_names]
    grouped = df.groupby(obs_keys)

    # Calculate mean log-normalized expression for each group.
    averaged = grouped[lognorm_cols].mean()

    # Calculate the percentage of cells in each group with expression above the
    # global median.
    # This is more robust to outliers than comparing means.
    median = df[lognorm_cols].median()
    percent_above = grouped[lognorm_cols].apply(lambda g: (g > median).mean() * 100)
    percent_above.columns = [v.replace("_lognorm", "") for v in percent_above.columns]

    # --- 4. Plotting ---
    if plot:
        # Create a single column for grouping to simplify plotting.
        group_col = "_group"
        if len(obs_keys) > 1:
            # If multiple keys, concatenate them into a single string.
            df[group_col] = df[obs_keys].astype(str).agg("_".join, axis=1)
        else:
            group_col = obs_keys[0]

        # "Melt" the DataFrame from wide to long format, which is required by seaborn.
        melted = df.melt(
            id_vars=[group_col],
            value_vars=lognorm_cols,
            var_name="gene",
            value_name="lognorm_expression",
        )
        melted["gene"] = melted["gene"].str.replace("_lognorm", "")

        logger.info("Generating plots…")
        for gene in var_names:
            plt.figure(figsize=(max(6, df[group_col].nunique() * 0.8), 6))
            data = melted[melted["gene"] == gene]
            plot_func = sns.violinplot if plot_type == "violin" else sns.boxplot
            plot_func(data=data, x=group_col, y="lognorm_expression", palette="viridis")
            # Overlay a stripplot to show individual data points.
            sns.stripplot(
                data=data,
                x=group_col,
                y="lognorm_expression",
                color="black",
                size=2.5,
                alpha=0.4,
                jitter=True,
            )
            plt.title(f"Log-Normalized Expression of {gene}")
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            plt.grid(axis="y", linestyle="--", alpha=0.6)
            if save_dir:
                import os
                os.makedirs(save_dir, exist_ok=True)
                filepath = os.path.join(save_dir, f"expression_{gene}.png")
                plt.savefig(filepath, bbox_inches='tight')
                plt.close()
            else:
                plt.show()            

    # --- 5. Return Results ---
    # Return a structured dictionary conforming to the ExpressionResults TypedDict.
    return {
        "individual_values": df.set_index(adata.obs.index),
        "averaged_values": averaged,
        "percent_above_median": percent_above,
    }


def preprocess_adata(
    filename: str,
    batch_name: str,
    apply_qc_filtering: bool = True,
    min_cells_per_gene: int = 3,
    min_genes_per_cell: int = 200,
    counts_quantile_low: float = 0.05,
    counts_quantile_high: float = 0.95,
    mt_pct_quantile_high: float = 0.95,
    leiden_resolution: float = 0.4,
) -> ad.AnnData:
    """
    Reads 10x H5 data and performs a standard scanpy preprocessing workflow.

    Args:
        filename: Path to the 10x H5 file.
        batch_name: An identifier for this data batch.
        ... and other QC parameters ...

    Returns:
        A fully processed AnnData object.
    """
    # --- 1. Data Loading and Initial Setup ---
    logger.info(f"Reading data from {filename} for batch '{batch_name}'…")
    if filename.name.endswith(".h5"):
        adata = sc.read_10x_h5(filename)
    elif filename.name.endswith(".h5ad"):
        adata = sc.read_h5ad(filename)
    else: 
        adata = sc.read_10x_mtx(filename)
        
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs["batch"] = batch_name

    # --- 2. Quality Control (QC) ---
    # Identify mitochondrial genes using a regular expression.
    adata.var["mt"] = adata.var_names.str.match(r"^(MT|mt)-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    if apply_qc_filtering:
        logger.info("Applying QC filtering…")
        # Filter out genes found in too few cells and cells with too few genes.
        sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
        sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)

        # Flag but do not remove cells with outlier QC metrics
        adata.obs["qc_issue"] = "OK"
        low = np.quantile(adata.obs.total_counts, counts_quantile_low)
        high = np.quantile(adata.obs.total_counts, counts_quantile_high)
        mt_high = np.quantile(adata.obs.pct_counts_mt, mt_pct_quantile_high)

        adata.obs.loc[adata.obs.total_counts < low, "qc_issue"] = "low_counts"
        adata.obs.loc[adata.obs.total_counts > high, "qc_issue"] = "high_counts"
        # Potential doublets
        adata.obs.loc[adata.obs.pct_counts_mt > mt_high, "qc_issue"] = "high_mt"
        # Stressed or dying cells

    # --- 3. Normalization and Data Storage ---
    # Store the original raw counts in a separate layer before normalization.
    adata.layers["counts"] = adata.X.copy()
    # Normalize by library size and apply log1p transformation.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Store the log-normalized data in another layer.
    adata.layers["l1p"] = adata.X.copy()
    # Store a snapshot of the full, processed data before any potential gene filtering.
    # This is crucial for functions like cell cycle scoring that need the full gene set.
    adata.raw = adata[:, :]

    # --- 4. Feature Scaling and Dimensionality Reduction ---
    # Scale data to unit variance and zero mean. Clip at max_value=10.
    sc.pp.scale(adata, max_value=10)
    logger.info("Running PCA, neighbors, UMAP, Leiden…")
    # Run the standard dimensionality reduction and clustering pipeline.
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_resolution)

    logger.info("Preprocessing complete.")
    return adata


def calculate_cell_cycle_scores(adata: ad.AnnData) -> ad.AnnData:
    """
    Calculates cell cycle phase scores based on Tirosh et al., 2015 gene lists.

    Args:
        adata: An AnnData object that has been preprocessed.

    Returns:
        The AnnData object, modified in-place with 'S_score', 'G2M_score', and 'phase'.
    """
    logger.info("Calculating cell cycle scores…")
    # This function requires the full gene set, so it must be run on an AnnData
    # object where the .raw attribute has been set.
    if adata.raw is None:
        raise ValueError("adata.raw is missing. Please run preprocessing first.")

    # Create a temporary AnnData object from the raw data. This avoids altering
    # the main adata object which might be subset to highly variable genes.
    adata_cc = adata.raw.to_adata()
    # Scaling is required for the scoring function.
    sc.pp.scale(adata_cc)

    # Filter the master gene lists to only include genes present in our data.
    s_genes = [g for g in S_GENES if g in adata_cc.var_names]
    g2m_genes = [g for g in G2M_GENES if g in adata_cc.var_names]

    if not s_genes or not g2m_genes:
        logger.warning("Not all cell cycle genes were found in the data.")

    # Use scanpy's built-in function to calculate S and G2/M phase scores.
    sc.tl.score_genes_cell_cycle(adata_cc, s_genes=s_genes, g2m_genes=g2m_genes)

    # Copy the calculated scores
    # and phase predictions back to the original AnnData object.
    adata.obs["S_score"] = adata_cc.obs["S_score"]
    adata.obs["G2M_score"] = adata_cc.obs["G2M_score"]
    adata.obs["phase"] = adata_cc.obs["phase"]

    logger.info("Cell cycle scoring complete.")
    return adata

def plot_composition(
    adata: ad.AnnData,
    groupby: Union[str, List[str]],
    split_by: Union[str, List[str]],
    normalize: bool = True,
    figsize: tuple = (10, 6),
    colormap: str = "viridis",
    save_path: Union[str, None] = None
):
    if isinstance(groupby, str):
        groupby = [groupby]
    if isinstance(split_by, str):
        split_by = [split_by]

    keys = groupby + split_by
    _ensure_in(adata.obs.columns, keys, "Group keys")

    df = adata.obs[keys].copy()
    df["count"] = 1
    composition = df.groupby(groupby + split_by).count().reset_index()
    composition = composition.pivot(index=groupby, columns=split_by, values="count").fillna(0)

    if normalize:
        composition = composition.div(composition.sum(axis=1), axis=0) * 100
        ylabel = "Percentage of Cells"
    else:
        ylabel = "Cell Count"

    composition.plot(kind="bar", stacked=True, figsize=figsize, colormap=colormap)
    plt.ylabel(ylabel)
    plt.xlabel("_".join(groupby))
    plt.xticks(rotation=45, ha="right")
    plt.title("Cell Composition by Group")
    plt.legend(title="_".join(split_by), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    else:
        plt.show()
