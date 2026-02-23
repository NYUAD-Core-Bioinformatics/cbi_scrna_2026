"""
sctoolkit: A collection of helper functions to speed up scanpy analysis.
"""

__version__ = "0.01"

# Import the main analysis functions from the helpers module to make them
# directly available at the package level.
# You can also expose the gene lists if you think they might be useful
# for users to inspect or use in their own custom functions.
from .helpers import (
    G2M_GENES,
    S_GENES,
    analyze_expression,
    calculate_cell_cycle_scores,
    preprocess_adata,
    plot_composition,
)
