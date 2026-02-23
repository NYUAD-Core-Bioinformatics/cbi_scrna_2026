# sctoolkit

`sctoolkit` is a Python package that provides a collection of helper functions designed to speed up and streamline common workflows in single-cell analysis with Scanpy.

## Key Features

* **Streamlined Preprocessing**: Quickly load 10x Genomics H5 data and run a standard preprocessing workflow with a single function.
* **Expression Analysis**: Easily analyze and visualize gene expression across different cell clusters or conditions.
* **Cell Cycle Scoring**: Calculate cell cycle scores using the established gene sets from Tirosh et al., 2015.

## Installation

You can install `sctoolkit` directly from the GitHub repository.

**Prerequisites:**

* Python 3.10+
* `pip` and `git`

**Install from GitHub**

```bash
# Create a new environment
conda create --name sctoolkit python=3.11 -y

# Activate it
conda activate sctoolkit

# Install scanpy and its core dependencies from conda-forge
conda install -c conda-forge scanpy python-igraph leidenalg jupyterlab -y

# Install any remaining dependencies using pip
pip install harmonypy scvi-tools seaborn

# Install sctoolkit
git clone https://github.com/gs512/das2_sc_toolkit.git
cd das2_sc_toolkit
pip install -e .

```

## License

This project is licensed under the MIT License.
