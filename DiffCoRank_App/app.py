import streamlit as st
import pandas as pd
import numpy as np
import os
import networkx as nx
import matplotlib.pyplot as plt
from io import BytesIO
import time
import shutil
import streamlit.components.v1 as components


from diffcorank_core import (
    columnFilter_loc, fetch_data_concurrently, update_extName,
    filter_data_by_iteration, generate_non_zero_mask,
    filter_data_by_moresample, filter_extName_empty_or_space,
    filter_and_compare_gene_df, compute_correlations,
    find_strongly_connected_genes, find_strongly_connected_genes_subset,
    compute_adjacency, quick_umap, density_cluster,
    analyze_module, write_compare_modules
)

st.set_page_config(
    page_title="DiffCoRank",
    page_icon="🧬",
    layout="wide"
)


st.title("DiffCoRank")
st.subheader("An Interactive Pipeline for Differential Co-expression Analysis")

st.sidebar.header("About")
st.sidebar.markdown("""
**Developed by:** Anirban Chakraborty<br>
**Contact:** chakra96@msu.edu<br>
**Version:** 1.0<br>
**All Rights Reserved © 2025**
""", unsafe_allow_html=True)

st.sidebar.link_button("GitHub Repo 🐱", "https://github.com/msureil/DiffCoRank")



@st.cache_data
def fetch_metadata(gene_ids):
    return fetch_data_concurrently(gene_ids)

for flag in ['filtered','correlation_done','scg_done','clustering_done']:
    if flag not in st.session_state:
        st.session_state[flag] = False


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


rlog   = pd.read_csv(rlog_file, sep=None, engine='python', index_col=0)
raw    = pd.read_csv(raw_file, sep=None, engine='python', index_col=0)
sample = pd.read_csv(sample_file)


name = raw.index.astype(str).tolist() #Row Names
nc1 = columnFilter_loc(rlog, sample, 'source', 'far', 'sample')  # Normalized Count for 'far'
nc2 = columnFilter_loc(rlog, sample, 'source', 'near', 'sample')  # Normalized Count for 'near'
rc1 = columnFilter_loc(raw, sample, 'source', 'far', 'sample')  # Raw Count for 'far'
rc2 = columnFilter_loc(raw, sample, 'source', 'near', 'sample')  # Raw Count for 'near'

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


with st.spinner("Fetching gene metadata..."):
    #gene_df = fetch_metadata(data['name'])
    
    gene_df = pd.read_csv('ensembl_gene_info_batch_filtered.csv')
st.success(f"Fetched metadata for {len(gene_df)} genes.")


meta_csv = gene_df.to_csv(index=False).encode('utf-8')
st.download_button(
    label="Download gene metadata", data=meta_csv,
    file_name="gene_metadata.csv", mime="text/csv"
)


st.subheader("Filter Stage")
st.write("Adjust filter parameters in the left and click 'Run Filtering' to proceed.")


st.sidebar.divider()

with st.sidebar.form("filter_form"):
    st.header("Gene Filter Parameters")
    min_raw_count        = st.slider(
        "Min Raw Counts threshold", 0, 300, 100
    )
    min_normalized_count = st.slider(
        "Min Normalized Counts threshold", 0, 10, 0
    )
    min_samp             = st.slider(
        "Min samples per gene", 1, 20, 3
    )
    min_len              = st.number_input(
        "Min gene length (bp)", 50, 1000, 100
    )
    cor_fdr_levels       = st.multiselect(
        "FDR levels to compute",
        [0.1, 0.5, 0.9],
        default=[0.1, 0.5]
    )

    run_filter = st.form_submit_button("Run Filtering")

if run_filter:  

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

  
    st.subheader("Connectivity Scatter Plot (50% vs 10% FDR)")
    fig3, ax3 = plt.subplots()
    ax3.scatter(sum2, sum1, s=1)
    ax3.set_xlabel("Number of correlations above 50% FDR")
    ax3.set_ylabel("Number of correlations above 10% FDR")
    st.pyplot(fig3)

   
    st.subheader("FDR Thresholds")


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


    if table_data:
        summary_df = pd.DataFrame(table_data)
        st.table(summary_df.set_index("FDR Level"))
    else:
        st.warning("FDR threshold data is not available.")

    st.subheader("Select the Minimum number of correlations for SCG")
  
    def sync_to_num():
        st.session_state.scg_num_input = st.session_state.scg_slider

    def sync_to_slider():
        st.session_state.scg_slider = st.session_state.scg_num_input


    if 'scg_slider' not in st.session_state:
        st.session_state.scg_slider           = 200
        st.session_state.scg_num_input        = 200
        st.session_state.sc_threshold_applied = False


    col1, col2 = st.columns([3,1])
    with col1:
        st.slider(
            "Min connections for SCG",
            1, 800,
            key="scg_slider",
            on_change=sync_to_num,
            label_visibility="collapsed"
        )
    with col2:
        st.number_input(
            "Value",
            1, 800,
            key="scg_num_input",
            on_change=sync_to_slider,
            label_visibility="collapsed"
        )


    if st.button("Apply Threshold"):
        st.session_state.sc_threshold     = st.session_state.scg_slider
        st.session_state.sc_threshold_applied = True


    if st.session_state.sc_threshold_applied:
        sc_threshold = st.session_state.sc_threshold
        st.success(f"The selected correlation threshold is: **{sc_threshold}** for 50% FDR")

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

    st.sidebar.divider()
    with st.sidebar.form("clustering_form"):
        st.header("Clustering Parameters")
        umap_n_neighbors = st.slider(
            "UMAP n_neighbors", 1, 20, 4
        )
        umap_min_dist = st.number_input(
            "UMAP min_dist", 0.00, 0.10, 0.005, 0.001, format="%.3f"
        )
        eps = st.number_input(
            "DBSCAN eps", 0.1, 1.0, 0.25, 0.1
        )
        min_samples = st.number_input(
            "DBSCAN min_samples", 1, 40, 14
        )
        min_mod_size = st.number_input(
            "Min module size", 1, 70, 50
        )
        run_clustering = st.form_submit_button("Run Clustering Pipeline")

    if run_clustering:
        st.session_state.clustering_done = True
        st.success("Clustering Pipeline Started Successfully 🎉")
    
if "download_done" not in st.session_state:
    st.session_state.download_done = False


@st.cache_data(show_spinner=False)
def get_adjacency(corrsSC, cap=1, threshold=0.03):
    return compute_adjacency(corrsSC, cap=cap, threshold=threshold)

if st.session_state.clustering_done:
    # initialize wizard step
    if "clust_step" not in st.session_state:
        st.session_state.clust_step = 0

    step = st.session_state.clust_step

    # --- STEP 0: Adjacency & TOM ---
    if step == 0:
        st.header("Step 1: Adjacency & TOM Distance")
        with st.spinner("Mapping connections. This is the most computationally intensive step, please wait!.…"):
            pbar = st.progress(0)
            for pct in range(0, 40):
                pbar.progress(pct+1); time.sleep(0.03)
            adSC = get_adjacency(st.session_state.corrsSC)
            st.session_state.adSC = adSC
            for pct in range(40, 100):
                pbar.progress(pct+1); time.sleep(0.005)
            pbar.empty()

            fig = plt.figure()
            plt.hist(
                np.log(adSC['Adj'][adSC['Adj']>0].flatten()),
                bins=100
            )
            plt.title("Histogram of log-transformed adjacency")
            st.pyplot(fig)

        if st.button("→ Next: UMAP"):
            st.session_state.clust_step += 1
            st.rerun()

    # --- STEP 1: UMAP ---
    elif step == 1:
        st.header("Step 2: UMAP of Strongly Connected Genes")
        umapSC = quick_umap(
            st.session_state.adSC['TOMdist'],
            n_neighbors=umap_n_neighbors,
            min_dist=umap_min_dist
        )
        st.session_state.umapSC = umapSC

        fig = plt.figure()
        plt.scatter(umapSC['x'], umapSC['y'], s=6)
        plt.title("UMAP (TOM Distance)")
        st.pyplot(fig)

        cols = st.columns([1,1,1])
        if cols[0].button("← Back"):
            st.session_state.clust_step -= 1; st.rerun()
        if cols[2].button("→ Next: Modules"):
            st.session_state.clust_step += 1; st.rerun()

    # --- STEP 2: DBSCAN Modules ---
    elif step == 2:
        st.header("Step 3: DBSCAN Module Detection")
        modsSC = density_cluster(
            st.session_state.umapSC,
            eps=eps,
            min_samples=min_samples,
            min_mod_size=min_mod_size
        )
        st.session_state.modsSC = modsSC

        fig, ax = plt.subplots()
        for m in np.unique(modsSC):
            pts = st.session_state.umapSC[modsSC == m]
            label = f"{m+1}" if m >= 0 else "Noise"
            ax.scatter(pts['x'], pts['y'], s=6, label=label)
        ax.legend(markerscale=2)
        st.pyplot(fig)

        cols = st.columns([1,1,1])
        if cols[0].button("← Back"):
            st.session_state.clust_step -= 1; st.rerun()
        if cols[2].button("→ Next: Hub Summary"):
            st.session_state.clust_step += 1; st.rerun()

    # --- STEP 3: Hub Summary & Download / Reset ---
    elif step == 3:
        


        with st.expander("🔍 View Step 1: Adjacency & TOM", expanded=False):
            adSC = st.session_state.adSC
            fig = plt.figure()
            plt.hist(
                np.log(adSC['Adj'][adSC['Adj'] > 0].flatten()),
                bins=100
            )
            plt.title("Histogram of log-transformed adjacency values")
            plt.xlabel("Log(adj)")
            plt.ylabel("Frequency")
            st.pyplot(fig)

        # 2️⃣ UMAP
        with st.expander("🔍 View Step 2: UMAP", expanded=False):
            umapSC = st.session_state.umapSC
            fig = plt.figure()
            plt.scatter(umapSC['x'], umapSC['y'], s=6)
            plt.title("UMAP of Strongly Connected Genes")
            plt.xlabel("UMAP1")
            plt.ylabel("UMAP2")
            st.pyplot(fig)

        # 3️⃣ Modules (DBSCAN)
        with st.expander("🔍 View Step 3: Module Clustering", expanded=False):
            modsSC = st.session_state.modsSC
            umapSC = st.session_state.umapSC
            fig, ax = plt.subplots()
            for m in np.unique(modsSC):
                pts = umapSC[modsSC == m]
                label = f"{m+1}" if m >= 0 else "Noise"
                ax.scatter(pts['x'], pts['y'], s=6, label=label)
            ax.legend(markerscale=2)
            ax.set_title("DBSCAN Clustering of SCGs")
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            st.pyplot(fig)
        st.header("Step 4: Hub Gene Summary")

        if st.session_state.get("download_done", False):
            st.success("Download complete! 🎉")



        if 'summary_df' not in st.session_state:
            with st.spinner("Analyzing modules…"):
                mod_dir = "resultsandplots/modules"
                os.makedirs(mod_dir, exist_ok=True)
                for m in np.unique(st.session_state.modsSC):
                    if m >= 0:
                        analyze_module(
                            st.session_state.modsSC,
                            m,
                            st.session_state.final_data,
                            st.session_state.adSC,
                            mod_dir
                        )
                st.session_state.summary_df = write_compare_modules(mod_dir)


        df = st.session_state.summary_df.copy()
        if 'identifier' in df:
            df['Module Number'] = df['identifier']+1
        st.dataframe(df[['Module Number','hubGene']], hide_index=True)

        csv_bytes = df.to_csv(index=False).encode()
        if st.download_button(
            "📥 Download Full Hub Gene List (CSV)",
            data=csv_bytes,
            file_name="hub_summary.csv",
            key="hub_download"
        ):
            st.session_state.download_done = True
            st.success("Download complete! 🎉")

 
        with st.expander("⚠️ Danger Zone", expanded=False):
            if st.button("🔄 Reset All", key="reset_final"):
                st.session_state.clear()
                shutil.rmtree("resultsandplots/modules", ignore_errors=True)
                st.rerun()


