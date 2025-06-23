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
# 🐳 DiffCoRank Docker Deployment Guide

## Quick Start with Docker

Follow these steps to get DiffCoRank up and running using Docker.

1. **Pull the official image from Docker Hub**  
   ```bash
   docker pull anirban1231/diffcorank:latest
   ```

2. **Run the container**
   ```bash
   docker run -d -p 8501:8501 \
     --name diffcorank_app \
     anirban1231/diffcorank:latest
   ```

3. **Open in your browser**
   Navigate to:
   ```
   http://localhost:8501
   ```

4. **Stop & remove**
   ```bash
   docker stop diffcorank_app
   docker rm diffcorank_app
   ```

---

## 💻 Platform-Specific Instructions

### 🐧 Linux (Ubuntu/Debian)
```bash
# Install Docker
sudo apt update
sudo apt install -y docker.io
sudo systemctl enable --now docker

# (Optional) Allow your user to run Docker without sudo
sudo usermod -aG docker $USER
# Log out & back in, or run:
newgrp docker

# Pull & run DiffCoRank
docker pull anirban1231/diffcorank:latest
docker run -d -p 8501:8501 \
  --name diffcorank_app \
  anirban1231/diffcorank:latest
```

### 🍎 macOS
1. Download and install **Docker Desktop** from
   [https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop)
2. Open **Terminal** and run:
   ```bash
   docker pull anirban1231/diffcorank:latest
   docker run -d -p 8501:8501 \
     --name diffcorank_app \
     anirban1231/diffcorank:latest
   ```
3. Browse to:
   ```
   http://localhost:8501
   ```

### ⊞ Windows
#### Docker Desktop (Recommended)
1. Install **Docker Desktop** from
   [https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop)
2. Open **PowerShell** as Administrator and run:
   ```powershell
   docker pull anirban1231/diffcorank:latest
   docker run -d -p 8501:8501 `
     --name diffcorank_app `
     anirban1231/diffcorank:latest
   ```
3. Visit:
   ```
   http://localhost:8501
   ```

#### WSL2 Integration
1. Ensure **WSL2** is installed and Docker Desktop is configured to use the WSL2 engine
2. In your WSL2 distribution shell:
   ```bash
   docker pull anirban1231/diffcorank:latest
   docker run -d -p 8501:8501 \
     --name diffcorank_app \
     anirban1231/diffcorank:latest
   ```
3. Open `http://localhost:8501` in your Windows browser

---

## ⚙️ Local Build & Development (Optional)
```bash
git clone https://github.com/msureil/DiffCoRank.git
cd DiffCoRank/DiffCoRank_App

# Build & tag
docker build -t anirban1231/diffcorank:latest .

# Run
docker run -d -p 8501:8501 \
  --name diffcorank_app \
  anirban1231/diffcorank:latest
```

## 📚 User Guide

Follow these steps to load your data, run the full DiffCoRank pipeline, and see your results.

---

### 1. Prepare Your Input Files

- **Rlog Normalized Counts (CSV/TXT)**  
  A matrix of log-transformed, normalized expression values.  
  - **Rows:** Gene IDs  
  - **Columns:** Sample names  
- **Raw Counts (CSV/TXT)**  
  Your original raw gene-count matrix.  
  - Same shape and sample order as the Rlog file.  
- **Sample Metadata (CSV)**  
  Metadata linking each sample to an experimental condition.  
  - Must contain at least two columns:  
    1. `sample` (matching the column names in your counts files)  
    2. `source` (e.g. "near" / "far")

---

### 2. Launch the App

```bash
# If you're using Docker:
docker pull anirban1231/diffcorank:latest
docker run -d -p 8501:8501 --name diffcorank_app anirban1231/diffcorank:latest
open http://localhost:8501
```

Or run locally:

```bash
pip install -r requirements.txt
streamlit run app.py
```

---

### 3. Upload Your Data

1. In the **Upload Your Data** section:
   * Drag & drop (or click) to upload:
     * **Rlog Normalized Counts**
     * **Raw Counts**
     * **Sample Metadata**
2. Once all three files appear, the **Run Filtering** button will activate.

---

### 4. Filter Stage

1. In the sidebar under **Gene Filter Parameters**, adjust:
   * **Min Raw Counts threshold**
   * **Min Normalized Counts threshold**
   * **Min samples per gene**
   * **Min gene length (bp)**
2. Click **Run Filtering**.
3. Review the **Filtering Summary** table showing the number of genes kept/excluded.

---

### 5. Correlation & SCG Selection

1. View **Significant Correlations Distribution** histograms at your chosen FDR levels.
2. In **Connectivity Scatter Plot**, drag to inspect the spread.
3. Under **Select the Minimum number of correlations for SCG**:
   * Drag the slider or type a number.
   * Click **Apply Threshold** to lock it in.
   * Click **Run Strongly Connected Genes Selection**.
4. You'll see the count and percentage of SCGs.

---

### 6. Clustering Wizard

Use the **Next** and **Back** buttons to step through:

* **Step 1: Adjacency & TOM**
  The progress bar displays adjacency and TOM computation.
* **Step 2: UMAP**
  Visualise SCGs in 2D.
* **Step 3: DBSCAN Modules**
  Identify colored modules with DBSCAN.
* **Step 4: Hub Gene Summary**
  Review and export your hub-gene list.

---

### 7. Export & Reset

* In **Step 4**, click **📥 Download Full Hub Gene List (CSV)** to save your results.
* When you're done, open the **⚠️ Danger Zone** expander and click **🔄 Reset All** to clear state and start over.

---

### 🎥 Video Tutorial

Watch the full walkthrough here:
👉 [DiffCoRank Tutorial](https://youtu.be/YOUR_VIDEO_ID)
---


## 📖 Additional Resources

- **Requirements**: Python 3.8+, Streamlit, pandas, numpy, networkx, matplotlib
- **Configuration**: Customize theme or port in `.streamlit/config.toml`
- **Troubleshooting**:
  - Permission errors on Linux? Add user to `docker` group
  - Windows issues? Ensure virtualization is enabled in BIOS
  - Port conflict? Change `8501` to another port in both commands

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


