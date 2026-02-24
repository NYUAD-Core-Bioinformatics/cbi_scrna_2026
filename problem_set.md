# Problem Set — Robust pre-processing & differential expression (scRNA-seq)

This problem set follows the workflow from the slide deck (QC → Normalization → Feature selection → Dimensionality reduction → Neighbors graph → Clustering → Annotation → DGE). 

## Part A — QC thresholds (cell filtering)

The slide deck highlights four QC axes: mitochondrial reads, read counts per cell, genes detected per cell, and doublets.

1) **Mitochondrial fraction threshold**
   - (a) Make the standard QC plots (e.g., `total_counts` vs `n_genes_by_counts` colored by `pct_counts_mt`, and violin plots). 
   - (b) Explain *how* you chose your mitochondrial threshold (e.g., a percentile cutoff, an elbow in the distribution, or a robust outlier rule).
   - (c) Report the **exact numeric threshold** you used (e.g., `pct_counts_mt < 12.5`).

2) **Library size (total counts) thresholds**
   - (a) Describe the evidence you used from the distribution/scatterplots to identify low-count and high-count outliers.
   - (b) State whether you used **percentiles**. 
   - (c) Report the **exact numeric low and high cutoffs** you used for `total_counts`.

3) **Genes detected thresholds**
   - (a) Explain how you chose cutoffs for `n_genes_by_counts` (low-complexity vs potential doublets/high complexity).
   - (b) Report the **exact numeric cutoffs** you used.

4) **Doublet detection**
   - (a) Run a doublet detection method (e.g., Scrublet) and show how doublet scores distribute across the embedding.
   - (b) Explain how you turned doublet scores into a binary “doublet/not doublet” decision (threshold rule).
   - (c) Report the **exact threshold** used and the **fraction of cells removed**.
   - (d) Briefly describe how doublets can impact downstream results if left in. 

5) **QC decision summary**
   - Provide a small summary table with your chosen thresholds (one row per metric) and a 1–2 sentence justification per row.
   - Report how many cells remain after QC filtering.

## Part B — Normalization choices

The deck motivates normalizing by count depth and log-transforming. 

6) **Normalization**
   - (a) State the normalization method you used (e.g., `normalize_total` to a target sum).
   - (b) Explain what the normalization is correcting for (in plain language).
   - (c) Show one diagnostic figure comparing the distribution of expression values *before vs after* log transform.

## Part C — Clustering & embedding

7) **Highly-variable genes (HVGs)**
   - (a) How many HVGs did you select, and why?
   - (b) Provide a plot or description demonstrating the mean–variance relationship motivation for HVGs. 

8) **Dimensionality reduction and neighbors graph**
   - (a) How many PCs did you keep? Describe your rationale (e.g., elbow plot, variance explained, stability).
   - (b) In 2–3 sentences, describe what the neighbors graph is used for in the pipeline. 

9) **Clustering**
   - (a) What clustering method and resolution did you use (e.g., Leiden resolution)?
   - (b) Show a UMAP (or equivalent) labeled by clusters.
   - (c) Describe one way you checked whether clusters were biologically plausible (markers, composition, etc.).

## Part D — Differential expression (DE)

10) **Top DE genes per cluster**
   - Run DE to find marker genes for each cluster (either **Leiden clusters** or **annotated cell types**, but be consistent).
   - For **each cluster**, report the **top 10 DE genes** (as a table or clearly formatted list).
   - Briefly interpret 1–2 clusters: do the top markers match expected biology?

11) **Explain the DE results dataframe**
   - Using the dataframe returned by `sc.get.rank_genes_groups_df(...)`, describe **each column** (what it means and how it is computed/used).
   - At minimum, cover: `names`, `scores`, `logfoldchanges`, `pvals`, `pvals_adj`, and any `% expressed` columns if present.

12) **QC and DE**
   - In 2–4 sentences, explain why including low-QC cells can blur differential expression signals.