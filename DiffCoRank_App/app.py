import streamlit as st
import pandas as pd
import numpy as np
import os
import networkx as nx
import matplotlib.pyplot as plt
from io import BytesIO

# Import core functions
from diffcorank_core import (
    columnFilter_loc, fetch_data_concurrently, update_extName,
    filter_data_by_iteration, generate_non_zero_mask,
    filter_data_by_moresample, filter_extName_empty_or_space,
    filter_and_compare_gene_df, compute_correlations,
    find_strongly_connected_genes, find_strongly_connected_genes_subset,
    compute_adjacency, quick_umap, density_cluster,
    analyze_module, write_compare_modules
)

# Cache metadata fetch to avoid re-fetching
@st.cache_data
def fetch_metadata(gene_ids):
    return fetch_data_concurrently(gene_ids)

# Initialize session state flags
for flag in ['filtered','correlation_done','scg_done','clustering_done']:
    if flag not in st.session_state:
        st.session_state[flag] = False


# Page config
st.set_page_config(
    page_title="DiffCoRank",
    page_icon="🧬",
    layout="wide"
)

# Main title of the app
st.title("DiffCoRank")
st.subheader("An Interactive Pipeline for Differential Co-expression Analysis")

#About Section
st.sidebar.header("About")
st.sidebar.markdown("""
**Developed by:** Anirban Chakraborty<br>
**Contact:** chakra96@msu.edu<br>
**Version:** 1.0<br>
**All Rights Reserved © 2025**
""", unsafe_allow_html=True)

st.sidebar.link_button("GitHub Repo 🐱", "https://github.com/msureil/DiffCoRank")


with st.container(border=True):
    st.subheader("Upload Your Data")
    st.markdown("Begin the analysis by providing the three required files. Drag and drop or click to browse.")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("##### Normalized Counts")
        rlog_file = st.file_uploader(
            "Rlog Normalized Counts", 
            type=["csv", "txt"], 
            key="rlog",
            label_visibility="collapsed"
        )
        st.caption("File containing log-transformed, normalized expression data.")

    with col2:
        st.markdown("##### Raw Counts")
        raw_file = st.file_uploader(
            "Raw Expression Counts", 
            type=["csv", "txt"], 
            key="raw",
            label_visibility="collapsed"
        )
        st.caption("Raw gene count data from your experiment.")

    with col3:
        st.markdown("##### Sample Information")
        sample_file = st.file_uploader(
            "Sample Metadata", 
            type=["csv"], 
            key="sample",
            label_visibility="collapsed"
        )
        st.caption("File linking sample IDs to experimental conditions (e.g.,'near','far').")


if not (rlog_file and raw_file and sample_file):
    st.info("Please upload all three files to proceed with the analysis.")
    st.stop()

# Load data only once
rlog   = pd.read_csv(rlog_file, sep=None, engine='python', index_col=0)
raw    = pd.read_csv(raw_file, sep=None, engine='python', index_col=0)
sample = pd.read_csv(sample_file)


name = raw.index.astype(str).tolist() #Row Names
nc1 = columnFilter_loc(rlog, sample, 'source', 'far', 'sample')  # Normalized Count for 'far'
nc2 = columnFilter_loc(rlog, sample, 'source', 'near', 'sample')  # Normalized Count for 'near'
rc1 = columnFilter_loc(raw, sample, 'source', 'far', 'sample')  # Raw Count for 'far'
rc2 = columnFilter_loc(raw, sample, 'source', 'near', 'sample')  # Raw Count for 'near'

#Stored in data dictionary.
data = {
    "name": name,
    "extName": name,
    "nc1": nc1,
    "rc1": rc1,
    "nc2": nc2,
    "rc2": rc2
}

st.success("Data loaded and partitioned successfully.")
st.write(f"**Initial gene count:** {len(data['name'])}")

# --- Metadata Fetch and Download ---
with st.spinner("Fetching gene metadata..."):
    #gene_df = fetch_metadata(data['name'])
    gene_df = pd.read_csv('ensembl_gene_info_batch_filtered.csv')
st.success(f"Fetched metadata for {len(gene_df)} genes.")

# Download metadata button
meta_csv = gene_df.to_csv(index=False).encode('utf-8')
st.download_button(
    label="Download gene metadata", data=meta_csv,
    file_name="gene_metadata.csv", mime="text/csv"
)

# --- Filter Parameters ---
st.subheader("Filter Stage")
st.write("Adjust filter parameters in the left and click 'Run Filtering' to proceed.")

st.sidebar.divider() 
st.sidebar.header("Gene Filter Parameters")
min_raw_count = st.sidebar.slider("Min Raw Counts threshold", 0, 300, 100)
min_normalized_count = st.sidebar.slider("Min Normalized Counts threshold", 0, 10, 0)
min_samp      = st.sidebar.slider("Min samples per gene", 1, 20, 3)
min_len       = st.sidebar.number_input("Min gene length (bp)", 50, 1000, 100)
cor_fdr_levels= st.sidebar.multiselect("FDR levels to compute", [0.1,0.5,0.9], default=[0.1,0.5,0.9])

if st.sidebar.button("Run Filtering"):  
    # Filtering pipeline
    updated_data = update_extName(data, gene_df)
    updated_data_new = {
        "name": gene_df['id'].tolist(),
        "extName": updated_data['extName'],
        "nc1": nc1,
        "rc1": rc1,
        "nc2": nc2,
        "rc2": rc2
    }

    lengths = {
        'name': len(updated_data_new['name']),
        'extName': len(updated_data_new['extName']),
        'rc1': len(updated_data_new['rc1']),
        'rc2': len(updated_data_new['rc2']),
        'nc1': len(updated_data_new['nc1']),
        'nc2': len(updated_data_new['nc2']),
    }

    # Checking if all lengths are equal
    if len(set(lengths.values())) != 1:
        st.warning("Mismatch detected in data field lengths", icon="⚠️")


    updated_data_new['rc1'] = updated_data_new['rc1'].apply(pd.to_numeric, errors='coerce')
    updated_data_new['rc2'] = updated_data_new['rc2'].apply(pd.to_numeric, errors='coerce')
    updated_data_new['nc1'] = updated_data_new['nc1'].apply(pd.to_numeric, errors='coerce')
    updated_data_new['nc2'] = updated_data_new['nc2'].apply(pd.to_numeric, errors='coerce')

    total_raw_counts = updated_data_new['rc1'].sum(axis=1, skipna=True) + updated_data_new['rc2'].sum(axis=1, skipna=True)
    total_normalized_counts = updated_data_new['nc1'].sum(axis=1, skipna=True) + updated_data_new['nc2'].sum(axis=1, skipna=True)
    high_raw_counts = (total_raw_counts >= min_raw_count) & (total_normalized_counts > min_normalized_count)

    filtered_data_1count = filter_data_by_iteration(updated_data_new, high_raw_counts)
    data_f1 = filtered_data_1count ###

    s1 = (updated_data_new['rc1'] > 0).sum(axis=1)
    s2 = (updated_data_new['rc2'] > 0).sum(axis=1)
    more_samples = (s1 >= min_samp) & (s2 >= min_samp)
    fewsample_mask = generate_non_zero_mask(filtered_data_1count, filtered_data_1count['name'], min_samp)
    mask    = fewsample_mask ###
    filtered_data_2moresample = filter_data_by_moresample(filtered_data_1count, fewsample_mask)
    data_f2 = filtered_data_2moresample ###
    filtered_data_3noname = filter_extName_empty_or_space(filtered_data_2moresample)

    data_f3 = filtered_data_3noname ###
    filtered_data_4_otherfilters, filtered_gene_names_final = filter_and_compare_gene_df(filtered_data_3noname, gene_df, min_len)
    final_data = filtered_data_4_otherfilters
    st.session_state.final_data = filtered_data_4_otherfilters
    st.session_state.filt_genes = filtered_gene_names_final
    st.session_state.filtered = True
    initial_count = len(updated_data_new['name'])
    st.session_state.initial_count = initial_count
    filter1_count = len(filtered_data_1count['name'])
    st.session_state.filter1_count = filter1_count 
    filter2_count = len(filtered_data_2moresample['name'])
    st.session_state.filter2_count = filter2_count
    filter3_count = len(filtered_data_3noname['name'])
    st.session_state.filter3_count = filter3_count
    filter4_count = len(filtered_data_4_otherfilters['name'])
    st.session_state.filter4_count = filter4_count
    final_count = st.session_state.final_data['name'] 

# After filtering
if st.session_state.filtered:
    
    if st.session_state.get("filtered", True):
    
        initial_count = st.session_state.get('initial_count')
        filter1_count = st.session_state.get('filter1_count')
        filter2_count = st.session_state.get('filter2_count')
        filter3_count = st.session_state.get('filter3_count')
        final_count = len(st.session_state.final_data['name'])

        st.success(f"Filtering complete! Final gene count: {final_count}")

        st.subheader("Filtering Summary")

        tracking_log = {
            "Filter Stage": [
                "Initial Input",
                "After Raw/Normalized Count Threshold",
                "After Min Sample Count per Condition",
                "After Valid extName Check",
                "After Ensembl Annotation Filter"
            ],
            "Genes Remaining": [
                initial_count,
                filter1_count,
                filter2_count,
                filter3_count,
                final_count
            ]
        }

    
        summary_df = pd.DataFrame(tracking_log)

        genes_removed = [
            initial_count - filter1_count,
            filter1_count - filter2_count,
            filter2_count - filter3_count,
            filter3_count - final_count
        ]
        
        genes_removed.insert(0, np.nan)
        
        summary_df["Genes Removed in Step"] = genes_removed
        styled_df = summary_df.style.format({'Genes Remaining': '{:,}','Genes Removed in Step': lambda x: "—" if pd.isna(x) else f'{x:,.0f}'}).set_properties(**{'text-align': 'right'}, subset=['Genes Remaining', 'Genes Removed in Step'])
        st.table(styled_df)
        
        total_genes_excluded = initial_count - final_count
        st.info(f"Total genes remaining: **{final_count}**")
        st.warning(f"Total genes excluded: **{total_genes_excluded}**")

    # --- Correlation Stage ---
    if not st.session_state.correlation_done:
        with st.spinner("Computing Spearman correlations & FDR thresholds..."):
            results = compute_correlations(st.session_state.final_data)
            st.session_state.results = results
            st.session_state.correlation_done = True
    else:
        results = st.session_state.results

        
    
        # --- Correlation Exploration ---
    # Plot histograms of significant correlations per gene
    st.subheader("Significant Correlations Distribution")
    cstar10 = results['FDR_thresholds'][0.1]['cstar']
    cstar50 = results['FDR_thresholds'][0.5]['cstar']
    sum1 = ((results['C1_corr'] > cstar10).sum(axis=1) + (results['C2_corr'] > cstar10).sum(axis=1))
    sum2 = ((results['C1_corr'] > cstar50).sum(axis=1) + (results['C2_corr'] > cstar50).sum(axis=1))
    col1, col2 = st.columns(2)

    with col1:
        fig1, ax1 = plt.subplots()
        ax1.hist(sum1, bins=100, color='skyblue', edgecolor='black')
        ax1.set_title("Per Gene at 10% FDR")
        ax1.set_xlabel("Number of Significant Correlations")
        ax1.set_ylabel("Number of Genes")
        st.pyplot(fig1)

    with col2:
        fig2, ax2 = plt.subplots()
        ax2.hist(sum2, bins=100, color='salmon', edgecolor='black')
        ax2.set_title("Per Gene at 50% FDR")
        ax2.set_xlabel("Number of Significant Correlations")

        st.pyplot(fig2)

    # Scatter plot comparing connectivity
    st.subheader("Connectivity Scatter Plot (50% vs 10% FDR)")
    fig3, ax3 = plt.subplots()
    ax3.scatter(sum2, sum1, s=1)
    ax3.set_xlabel("Number of correlations above 50% FDR")
    ax3.set_ylabel("Number of correlations above 10% FDR")
    st.pyplot(fig3)

    
    st.subheader("FDR Thresholds")

    # Prepare data for the table
    table_data = []
    for lvl in cor_fdr_levels:
        thr = results['FDR_thresholds'].get(lvl, {})
        p_value = thr.get('pstar', 0)
        corr_value = thr.get('cstar', 0)
        
        table_data.append({
            "FDR Level": f"{lvl:.0%}", 
            "P-value (p*)": f"{p_value:.5f}", 
            "Correlation (c*)": f"{corr_value:.3f}" 
        })

    # Create a DataFrame and display it as a static table
    if table_data:
        summary_df = pd.DataFrame(table_data)
        st.table(summary_df.set_index("FDR Level"))
    else:
        st.warning("FDR threshold data is not available.")

    # SCG threshold selection
    st.subheader("Select the Minimum number of correlations for SCG")

    def sync_number_input():
        """When the slider changes, this function updates the number input's value."""
        st.session_state.scg_num_input = st.session_state.scg_slider

    def sync_slider():
        """When the number input changes, this function updates the slider's value."""
        st.session_state.scg_slider = st.session_state.scg_num_input


    if 'scg_slider' not in st.session_state:
        st.session_state.scg_slider = 200
        st.session_state.scg_num_input = 200


    col1, col2 = st.columns([3, 1])

  
    with col1:
        st.slider(
            "Min connections for SCG",
            min_value=1,
            max_value=800,
            key="scg_slider",        
            on_change=sync_number_input, 
            label_visibility="collapsed"
        )

    with col2:
        st.number_input(
            "Value",
            min_value=1,
            max_value=800,
            key="scg_num_input",     
            on_change=sync_slider,   
            label_visibility="collapsed"
        )

    
    sc_threshold = st.session_state.scg_slider
    st.write(f"The selected correlation threshold is: **{sc_threshold}** for 50% FDR")

    if st.button("Run Strongly Connected Genes Selection"):
        corrs = {
        'C1': results['C1_corr'], 'C2': results['C2_corr']}

        testSC = find_strongly_connected_genes(
            corrs,
            cstar10,
            cstar50,
            sc_threshold
        )
        corrsSC = find_strongly_connected_genes_subset(corrs, testSC)
        st.session_state.testSC = testSC
        st.session_state.corrsSC = corrsSC
        st.session_state.scg_done = True

# After SCG selection
if st.session_state.scg_done:
    st.success(
    f"Strongly connected genes: {st.session_state.testSC.sum():,} of {len(st.session_state.testSC):,} "
    f"({(st.session_state.testSC.sum() / len(st.session_state.testSC) * 100) if len(st.session_state.testSC) > 0 else 0:.1f}%)"
    )

    st.subheader("Now Proceed to Clustering Stage")
    st.write("Adjust Clustering parameters in the left and click 'Run Clustering' to proceed.")
    
    # Sidebar: clustering parameters
    st.sidebar.divider()
    st.sidebar.header("Clustering Parameters")
    umap_n_neighbors = st.sidebar.slider("UMAP n_neighbors", 1, 20, 4)
    umap_min_dist    = st.sidebar.number_input("UMAP min_dist", 0.00, 0.10, 0.005,0.001,format="%.3f")
    eps              = st.sidebar.number_input("DBSCAN eps", 0.1, 1.0, 0.25, 0.1)
    min_samples      = st.sidebar.number_input("DBSCAN min_samples", 1, 40, 14)
    min_mod_size     = st.sidebar.number_input("Min module size", 1, 70, 50)
    if st.sidebar.button("Run Clustering Pipeline"):
        st.session_state.clustering_done = True
    
# Clustering stage
if st.session_state.clustering_done:
    tabs = st.tabs(["Adjacency & TOM","UMAP","Modules","Hub Summary"])

    # Adjacency & TOM
    with tabs[0]:
        with st.spinner("Computing Adjacency and TOM Distance (This will take some moment. Sit tight!)..."):
            adSC = compute_adjacency(st.session_state.corrsSC, cap=1, threshold=0.03)
            st.session_state.adSC = adSC
            fig = plt.figure()
            plt.hist(np.log(adSC['Adj'][adSC['Adj']>0].flatten()), bins=100)
            plt.title("Histogram of log-transformed adjacency values")
            plt.xlabel("Log(adj)")
            plt.ylabel("Frequency")
            st.pyplot(fig)

    # UMAP
    with tabs[1]:
        with st.spinner("Running UMAP..."):
            umapSC = quick_umap(st.session_state.adSC['TOMdist'], n_neighbors=umap_n_neighbors, min_dist=umap_min_dist)
            st.session_state.umapSC = umapSC
            fig = plt.figure()
            plt.scatter(umapSC['x'], umapSC['y'], s=6)
            plt.title("UMAP of Strongly Connected Genes (TOM Distance)")
            plt.xlabel("UMAP1")
            plt.ylabel("UMAP2")
            st.pyplot(fig)

    # In your "Modules" tab
    with tabs[2]:
        with st.spinner("Running DBSCAN to find modules..."):
            modsSC = density_cluster(st.session_state.umapSC, eps=eps, min_samples=min_samples, min_mod_size=min_mod_size)
            st.session_state.modsSC = modsSC
            fig, ax = plt.subplots()
            unique_mods = np.unique(modsSC)
            
            
            for mod_id in unique_mods:
                label = f"{mod_id + 1}" if mod_id != -1 else "Noise"
                subset = st.session_state.umapSC[modsSC == mod_id]
                ax.scatter(subset['x'], subset['y'], label=label, s=6)
                
            ax.set_title('DBSCAN Clustering of Strongly Connected Genes')
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.legend(markerscale=2)
            
            plt.tight_layout()
            st.pyplot(fig)

    # Hub Summary
    with tabs[3]:
        with st.spinner("Analyzing modules and identifying hub genes..."):
            mod_dir = "resultsandplots/modules"
            os.makedirs(mod_dir, exist_ok=True)
            for m in np.unique(st.session_state.modsSC):
                if m >= 0:
                    _ = analyze_module(st.session_state.modsSC, m, st.session_state.final_data, st.session_state.adSC, mod_dir)
            
            summary_file = os.path.join(mod_dir, "hub_genes_list.csv")
            write_compare_modules(mod_dir, summary_file)

            summary_df = pd.read_csv(summary_file)

            if 'identifier' in summary_df.columns:

                summary_df['identifier'] = summary_df['identifier'] + 1
                
                
                summary_df.rename(columns={'identifier': 'Module Number'}, inplace=True)
            else:
                st.warning("Warning: 'identifier' column not found in summary file. Cannot rename or re-index.")

            columns_to_display = ['Module Number', 'hubGene']
            display_df = summary_df[columns_to_display]
            styled_df = display_df.style.set_properties(**{'text-align': 'center'})
            st.dataframe(styled_df, hide_index=True)
            full_mod_dir = os.path.abspath(mod_dir)
            st.success("✔ All module files and the List of the hub genes have been saved to the directory below:")
            st.code(full_mod_dir, language='bash')
