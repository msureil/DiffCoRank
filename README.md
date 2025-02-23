## DifCoRank: A Comprehensive Framework for Discovering Hub Genes and Differential Gene Co-expression in Brain Implant-Associated Tissue Responses

**Author**: Anirban Chakraborty<sup>1,3</sup> (chakra96@msu.edu), Erin K. Purcell<sup>1,2,3</sup>, Michael G. Moore<sup>3</sup>  
**Affiliations**:  
1. Department of Electrical and Computer Engineering, Michigan State University  
2. Department of Biomedical Engineering, Michigan State University  
3. Institute for Quantitative Health Science and Engineering, Michigan State University  

---

### Overview

**DifCoRank** is an integrated framework to identify differentially coexpressed gene modules and prioritize key regulators (hub genes) associated with tissue responses to implanted devices. The pipeline combines:

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

### License

This project is licensed under the [MIT License](LICENSE).

---

### Contact / Support

- **Anirban Chakraborty**: chakra96@msu.edu  
- **Erin K. Purcell**: epurcell@msu.edu  
- **Michael G. Moore**: moorem38@msu.edu

For questions, bug reports, or feature requests, please open an issue in this repository or email the authors.

---

### Acknowledgments

- This research was conducted at Michigan State University’s Department of Electrical & Computer Engineering, Department of Biomedical Engineering, and the Institute for Quantitative Health Science & Engineering.


