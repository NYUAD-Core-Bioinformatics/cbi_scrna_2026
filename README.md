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

# Install sctoolkit
git clone https://github.com/gs512/cbi_scrna_2026.git
cd cbi_scrna_2026
pip install -e .

```

**Datasets:**
[Download link](https://drive.google.com/drive/folders/12lpQN-rJMFK-hg0W1IMcDV8acZ4ul3dG?usp=sharing)

## License

This project is licensed under the MIT License.
