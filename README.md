## DifCoRank: A Comprehensive Framework for Discovering Hub Genes and Differential Gene Co-expression in Brain Implant-Associated Tissue Responses

**Author**: Anirban Chakraborty<sup>1,3</sup> (chakra96@msu.edu), Erin K. Purcell<sup>1,2,3</sup>, Michael G. Moore<sup>3</sup>  
**Affiliations**:  
1. Department of Electrical and Computer Engineering, Michigan State University  
2. Department of Biomedical Engineering, Michigan State University  
3. Institute for Quantitative Health Science and Engineering, Michigan State University  

---

### Overview

**DifCoRank** is an integrated framework designed to identify differentially coexpressed gene modules and prioritize key regulators (hub genes) associated with tissue responses to implanted devices. The pipeline combines:

- **RNA-Seq Data Preprocessing & Gene Filtering**  
- **Correlation-Based Module Identification**  
- **Hybrid Clustering** using UMAP + DBSCAN  
- **Multi-criteria Hub Gene Ranking** with centrality metrics (degree, closeness, betweenness, eigenvector)

This repository accompanies our paper:  
> Chakraborty, A., Purcell, E.K., & Moore, M.G. *DifCoRank: A Comprehensive Framework for Discovering Hub Genes and Differential Gene Co-expression in Brain Implant-Associated Tissue Responses*. (Under Review in BMC Bioinformatics, 2025).

---

### Features

- **Gene-Centric Approach**: Focuses on individual gene-gene correlations rather than only module-level analysis.  
- **UMAP + DBSCAN**: Enhanced clustering to handle non-linear high-dimensional data, revealing subtle co-expression structures.  
- **Network Centrality-Based Hub Detection**: Ranks genes by multiple centrality criteria, revealing key regulators of tissue responses.

---

### Contents

- `src/`: Source code for the DifCoRank pipeline.  
- `data/`: Example datasets (or placeholders) demonstrating input format.  
- `notebooks/`: Jupyter notebooks showing step-by-step usage and visualization examples.  
- `results/`: Example outputs (hub gene lists, cluster assignments, figures).  
- `LICENSE`: License details.  
- `README.md`: This file.

---

### Installation

1. **Clone the Repository**  
   ```bash
   git clone https://github.com/msureil/DifCoRank.git
   cd DifCoRank
   ```

2. **Create/Activate a Python Environment (Recommended)**  
   ```bash
   python -m venv venv
   source venv/bin/activate  # or venv\Scripts\activate on Windows
   ```

3. **Install Dependencies**  
   ```bash
   pip install -r requirements.txt
   ```
   Ensure your `requirements.txt` includes libraries such as `pandas`, `numpy`, `scipy`, `networkx`, `scikit-learn`, and `umap-learn`.

---

### Usage

Below is a high-level workflow to run DifCoRank on your dataset. More detailed instructions can be found in `notebooks/Demo.ipynb`.

1. **Prepare Your RNA-Seq Data**  
   - Provide raw or normalized expression counts in a CSV or TSV file, with rows as genes and columns as samples.  

2. **Run the Pipeline**  

3. **Review Outputs**  
   - **Clustered Modules**: Gene assignments to each co-expression cluster.  
   - **Hub Gene Ranks**: CSV files listing genes by average centrality rank.  
   - **Figures**: Plots of UMAP embeddings, DBSCAN clusters, and adjacency histograms.


### Citation

### License

This project is licensed under the [MIT License](LICENSE).

---

### Contact / Support

- **Anirban Chakraborty**: chakra96@msu.edu  
- **Erin K. Purcell**: epurcell@msu.edu  
- **Michael G. Moore**: moorem38@msu.edu

For questions, bug reports, or feature requests, please open an issue in this repository.

---

### Acknowledgments

- This research was conducted at Michigan State University’s Department of Electrical & Computer Engineering, Department of Biomedical Engineering, and the Institute for Quantitative Health Science & Engineering.


