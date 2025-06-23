## DiffCoRank: A Comprehensive Framework for Discovering Hub Genes and Differential Gene Co-expression in Brain Implant-Associated Tissue Responses

**Author**: Anirban Chakraborty<sup>1,3</sup> (chakra96@msu.edu), Erin K. Purcell<sup>1,2,3</sup>, Michael G. Moore<sup>3</sup>  
**Affiliations**:  
1. Department of Electrical and Computer Engineering, Michigan State University  
2. Department of Biomedical Engineering, Michigan State University  
3. Institute for Quantitative Health Science and Engineering, Michigan State University  

---

### Overview

**DiffCoRank** is an integrated framework to identify differentially coexpressed gene modules and prioritize key regulators (hub genes) associated with tissue responses to implanted devices. The pipeline combines:

- **RNA-Seq Data Preprocessing & Gene Filtering**  
- **Correlation-Based Module Identification**  
- **Hybrid Clustering** using UMAP + DBSCAN  
- **Multi-criteria Hub Gene Ranking** with centrality metrics (degree, closeness, betweenness, eigenvector)

  ![Flowchart](https://github.com/user-attachments/assets/263d36ee-6eb8-488b-b527-ab2bc3b42ff1)


This repository accompanies our paper:  
> Chakraborty, A., Purcell, E.K., & Moore, M.G. *DiffCoRank: A Comprehensive Framework for Discovering Hub Genes and Differential Gene Co-expression in Brain Implant-Associated Tissue Responses*. (Under Review in BMC Bioinformatics, 2025).

---

### Features

- **Gene-Centric Approach**: Focuses on individual gene-gene correlations rather than only module-level analysis.  
- **UMAP + DBSCAN**: Enhanced clustering to handle non-linear high-dimensional data, revealing subtle co-expression structures.  
- **Network Centrality-Based Hub Detection**: Ranks genes by multiple centrality criteria, revealing key regulators of tissue responses.

---
## Quick Start with Docker

1. **Pull the official image from docker hub**  
   ```bash
   docker pull anirban1231/diffcorank:latest
2. **Run the container**
    ```bash
    docker run -d -p 8501:8501 --name diffcorank_app anirban1231/diffcorank:latest
3. **Open in your browser**
    Navigate to: http://localhost:8501

4. **Stop & remove**
    ```bash
   docker stop diffcorank_app
   docker rm  diffcorank_app
---
## 💻 Platform-Specific Instructions

1.  Linux (Ubuntu/Debian)
    ```bash 
    # Install Docker
    sudo apt update
    sudo apt install -y docker.io
    sudo systemctl enable --now docker

    # (Optional) Manage permissions
    sudo usermod -aG docker $USER
    # Log out & back in, or run: newgrp docker

    # Pull & run
    docker pull anirban1231/diffcorank:latest
    docker run -d -p 8501:8501 --name diffcorank_app anirban1231/diffcorank:latest

2.  macOS
    - Install Docker Desktop from https://www.docker.com/products/docker-desktop
    - Open Terminal:
    ```bash
    docker pull anirban1231/diffcorank:latest
    docker run -d -p 8501:8501 --name diffcorank_app anirban1231/diffcorank:latest

  -Browse to http://localhost:8501
4.  Windows
    - Install Docker Desktop from https://www.docker.com/products/docker-desktop
    - Open PowerShell (Admin):
    ```bash
    docker pull anirban1231/diffcorank:latest
    docker run -d -p 8501:8501 --name diffcorank_app anirban1231/diffcorank:latest
    - Visit http://localhost:8501

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


