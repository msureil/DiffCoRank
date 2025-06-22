# diffcorank_core.py
"""
Core library for DiffCoRank analysis pipeline.
Provides data filtering, correlation computation, module detection,
and network-based hub gene analysis without any UI dependencies.
"""

import os
import json
import pandas as pd
import numpy as np
import requests
from concurrent.futures import ThreadPoolExecutor
import networkx as nx
from scipy.stats import spearmanr
from itertools import product
import umap
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score, davies_bouldin_score
from scipy.spatial.distance import pdist, squareform
from typing import Tuple
from io import BytesIO
import glob

# ----------------------
# Utility Functions
# ----------------------

def column_filter_loc(df_to_filter: pd.DataFrame,
                      metadata: pd.DataFrame,
                      group_col: str,
                      group_val: str,
                      sample_col: str) -> pd.DataFrame:
    """
    Keep columns of df_to_filter where metadata[group_col] == group_val.
    """
    samples = metadata.loc[metadata[group_col] == group_val, sample_col]
    cols = df_to_filter.columns.intersection(samples)
    return df_to_filter.loc[:, cols]

# Backward compatibility alias
columnFilter_loc = column_filter_loc


def split_batches(gene_ids: list, batch_size: int):
    """Yield slices of gene_ids of length batch_size."""
    for i in range(0, len(gene_ids), batch_size):
        yield gene_ids[i:i + batch_size]


# ----------------------------
def fetch_gene_data_batch(batch: list,
                          url: str,
                          headers: dict) -> list:
    """
    Fetch metadata for a batch of gene IDs from Ensembl REST API.
    """
    data = {"ids": batch}
    resp = requests.post(url, json=data, headers=headers)
    results = []
    if resp.ok:
        # Ensure we have a dict, even if resp.json() is None
        info = resp.json() or {}
        for gid in batch:
            # Coerce missing entries to empty dict
            g = info.get(gid) or {}
            length = abs(g.get('end', 0) - g.get('start', 0))
            results.append({
                'id': gid,
                'hasEnsData': bool(g),
                'name': g.get('display_name', ''),
                'type': g.get('biotype', ''),
                'length': length
            })
    else:
        for gid in batch:
            results.append({
                'id': gid,
                'hasEnsData': False,
                'name': '',
                'type': '',
                'length': None
            })
    return results
# ----------------------------



def fetch_data_concurrently(gene_ids: list,
                            batch_size: int = 1000,
                            max_workers: int = 10,
                            url: str = 'https://rest.ensembl.org/lookup/id',
                            headers: dict = None) -> pd.DataFrame:
    """
    Retrieve gene metadata for all gene_ids in parallel.
    """
    if headers is None:
        headers = {'Content-Type': 'application/json'}
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(fetch_gene_data_batch, batch, url, headers)
                   for batch in split_batches(gene_ids, batch_size)]
    records = []
    for f in futures:
        records.extend(f.result())
    return pd.DataFrame(records)


def update_extName(data, gene_df):
    gene_map = dict(zip(gene_df['id'], gene_df['name']))
    for i, gene_id in enumerate(data['extName']):
        if gene_id in gene_map:
            data['extName'][i] = gene_map[gene_id]
    return data


def filter_data_by_iteration(data, high_raw_counts):
    filtered_name = []
    filtered_extName = []
    filtered_rc1 = []
    filtered_rc2 = []
    filtered_nc1 = []
    filtered_nc2 = []
    
    for i, gene_id in enumerate(data['name']):
        if high_raw_counts.get(gene_id, True):
            filtered_name.append(data['name'][i])
            filtered_extName.append(data['extName'][i])
            filtered_rc1.append(data['rc1'].loc[gene_id].values.tolist())
            filtered_rc2.append(data['rc2'].loc[gene_id].values.tolist())
            filtered_nc1.append(data['nc1'].loc[gene_id].values.tolist())
            filtered_nc2.append(data['nc2'].loc[gene_id].values.tolist())

    filtered_data_1count = {
        'name': filtered_name,
        'extName': filtered_extName,
        'rc1': filtered_rc1,
        'rc2': filtered_rc2,
        'nc1': filtered_nc1,
        'nc2': filtered_nc2
    }

    return filtered_data_1count


def generate_non_zero_mask(data, ids, min_samp):
    rc1_df = pd.DataFrame(data['rc1'])
    rc2_df = pd.DataFrame(data['rc2'])

    s1 = (rc1_df > 0).sum(axis=1)
    s2 = (rc2_df > 0).sum(axis=1)

    more_samples = (s1 >= min_samp) & (s2 >= min_samp)

    mask_series = pd.Series(more_samples.values, index=ids)
    mask_series.index.name = 'id'
    return mask_series


def filter_data_by_moresample(data, high_raw_counts):
    filtered_name = []
    filtered_extName = []
    filtered_rc1 = []
    filtered_rc2 = []
    filtered_nc1 = []
    filtered_nc2 = []

    for i, gene_id in enumerate(data['name']):  
        if high_raw_counts.loc[gene_id]: 
            filtered_name.append(gene_id)
            filtered_extName.append(data['extName'][i])
            filtered_rc1.append(data['rc1'][i])  
            filtered_rc2.append(data['rc2'][i])
            filtered_nc1.append(data['nc1'][i])
            filtered_nc2.append(data['nc2'][i])

    filtered_data = {
        'name': filtered_name,
        'extName': filtered_extName,
        'rc1': filtered_rc1,
        'rc2': filtered_rc2,
        'nc1': filtered_nc1,
        'nc2': filtered_nc2
    }
    return filtered_data


def filter_extName_empty_or_space(data):
    filtered_name = []
    filtered_extName = []
    filtered_rc1 = []
    filtered_rc2 = []
    filtered_nc1 = []
    filtered_nc2 = []
    for i, ext_name in enumerate(data['extName']):
        #if ext_name.strip():
        if not pd.isna(ext_name) and str(ext_name).strip():
            filtered_name.append(data['name'][i])
            filtered_extName.append(ext_name)
            filtered_rc1.append(data['rc1'][i])
            filtered_rc2.append(data['rc2'][i])
            filtered_nc1.append(data['nc1'][i])
            filtered_nc2.append(data['nc2'][i])
    filtered_data = {
        'name': filtered_name,
        'extName': filtered_extName,
        'rc1': filtered_rc1,
        'rc2': filtered_rc2,
        'nc1': filtered_nc1,
        'nc2': filtered_nc2
    }

    return filtered_data


def filter_and_compare_gene_df(data, gene_df, min_len):
    gene_df['length'] = pd.to_numeric(gene_df['length'], errors='coerce')
    gene_df_subset = gene_df[gene_df['name'].isin(data['extName'])]

    not_coding_gene = gene_df_subset['type'] != 'protein_coding'
    too_short = gene_df_subset['length'] < min_len
    no_data = ~gene_df_subset['hasEnsData']

    mask = no_data | not_coding_gene | too_short
    filtered_gene_names = gene_df_subset.loc[mask, 'name']

    filtered_indices = [i for i, ext_name in enumerate(data['extName']) 
                        if ext_name not in filtered_gene_names.values]

    filtered_data = {
        'name': [data['name'][i] for i in filtered_indices],
        'extName': [data['extName'][i] for i in filtered_indices],
        'rc1': [data['rc1'][i] for i in filtered_indices],
        'rc2': [data['rc2'][i] for i in filtered_indices],
        'nc1': [data['nc1'][i] for i in filtered_indices],
        'nc2': [data['nc2'][i] for i in filtered_indices]
    }

    
    excluded_gene_names = set(data['extName']).intersection(set(filtered_gene_names))
    print(f"Genes filtered in this step (non-coding, short, or missing data): {len(excluded_gene_names)}")

    return filtered_data, filtered_gene_names

# ----------------------
# Core Analysis Steps
# ----------------------


def compute_correlations(data):

    
    C1_corr, C1_pvals = spearmanr(data['nc1'], axis=1)
    C2_corr, C2_pvals = spearmanr(data['nc2'], axis=1)

    lower_indices = np.tril_indices(C1_corr.shape[0], k=-1)

    C1_pvals_lower = C1_pvals[lower_indices]
    C1_corr_lower = C1_corr[lower_indices]
    C2_pvals_lower = C2_pvals[lower_indices]
    C2_corr_lower = C2_corr[lower_indices]

    combined_pvals = np.concatenate((C1_pvals_lower, C2_pvals_lower))
    combined_cvals = np.concatenate((C1_corr_lower, C2_corr_lower))


    sorted_indices = np.argsort(combined_pvals)
    combined_pvals_sorted = combined_pvals[sorted_indices]
    combined_cvals_sorted = combined_cvals[sorted_indices]

    M = len(combined_pvals_sorted)
    FPR_levels = [0.1, 0.5, 0.9]
    fdr_results = {}


    for FPR in FPR_levels:
        kstar_indices = np.where(combined_pvals_sorted < (np.arange(1, M + 1) * FPR / M))[0]


        if len(kstar_indices) > 0:
            kstar = max(kstar_indices)
            prc = kstar / M
            pstar = combined_pvals_sorted[kstar]
            cstar = abs(combined_cvals_sorted[kstar])

            fdr_results[FPR] = {
                "kstar": kstar,
                "prc": prc,
                "pstar": pstar,
                "cstar": cstar
            }
        else:
            fdr_results[FPR] = {
                "kstar": None,
                "prc": None,
                "pstar": None,
                "cstar": None
            }

    return {
        "C1_corr": C1_corr,
        "C1_pvals": C1_pvals,
        "C2_corr": C2_corr,
        "C2_pvals": C2_pvals,
        "FDR_thresholds": fdr_results
    }


def find_strongly_connected_genes(corrs, cstar10, cstar50, min_corr_threshold):
    sum1 = np.sum(corrs['C1'] > cstar10, axis=1) + np.sum(corrs['C2'] > cstar10, axis=1)
    sum2 = np.sum(corrs['C1'] > cstar50, axis=1) + np.sum(corrs['C2'] > cstar50, axis=1)
    testSC = sum2 > min_corr_threshold
    return testSC

def find_strongly_connected_genes_subset(corrs, test):
    corrsSC = {
        'C1': corrs['C1'][np.ix_(test, test)],
        'C2': corrs['C2'][np.ix_(test, test)]
    }
    return corrsSC



def compute_TOM_distance(D):
    N = D.shape[0]
    np.fill_diagonal(D,0) 
    rowSum = np.sum(D,axis=1) 
    TOM = np.zeros((N, N), dtype=np.float64) 
    for i in range(N):
        for j in range(i+1, N):
            shared = 0.0
            for k in range(N):
                shared += D[i,k]*D[k,j]
            numerator = shared + D[i,j]     
            denominator = min(rowSum[i],rowSum[j]) + 1.0 - D[i,j] 
            if denominator!= 0:
                TOM[i,j] = numerator/denominator
            else:
                TOM[i,j] = 0.0
            TOM[j,i] = TOM[i,j]
    TOMdist = 1.0 - TOM   
    return TOMdist




def compute_adjacency(corrs, cap=1, threshold=0.03):

    C1 = corrs['C1'].copy()
    C2 = corrs['C2'].copy()
    C1[C1< threshold] = 0
    C2[C2< threshold] = 0
    A_diff = np.abs(C1 - C2)
    A = np.clip(A_diff,0,cap)/cap   
    np.fill_diagonal(A, 1.0)
    TOMdist = compute_TOM_distance(A)

    return {
        'cap': cap,
        'Adj': A,        #Adj Matrix
        'TOMdist': TOMdist  #NxN matrix of topological-overlap distances
    }


def quick_umap(distance_matrix, n_neighbors=15, min_dist=0.1):
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=2,
        random_state=123,   
        metric='precomputed'  
    )
    embedding = reducer.fit_transform(distance_matrix)

    
    return pd.DataFrame({'x': embedding[:, 0], 'y': embedding[:, 1]})
    


def density_cluster(umap_embedding, eps, min_samples, min_mod_size):
    
    if hasattr(umap_embedding, 'values'):
        coords = umap_embedding[['x', 'y']].values
    else:
        coords = umap_embedding

    #DBSCAN in Euclidean space. As I have used scikit-learn - it labels noise as -1
    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='euclidean')
    dbscan.fit(coords)
    mods = dbscan.labels_  #example: 0,1,2,... or -1 for noise

    unique_labels = set(mods) - {-1}    #(-1) is noise here. 
    for m in unique_labels:
        cluster_size = np.sum(mods == m)
        if cluster_size <= min_mod_size:
            mods[mods == m] = -1  #Mark these points as noise

    merge0=True
    if merge0:
        noise_indices = np.where(mods == -1)[0]
        if len(noise_indices) > 0:
            cluster_indices = np.where(mods != -1)[0]
            if len(cluster_indices) == 0:
                pass
            else:
                for idx in noise_indices:
                    distances = np.linalg.norm(coords[cluster_indices] - coords[idx], axis=1)
                    nearest_idx = cluster_indices[np.argmin(distances)]
                    mods[idx] = mods[nearest_idx]

    valid_labels = sorted(label for label in np.unique(mods) if label != -1)

    label_map = {}
    for i, old_label in enumerate(valid_labels):
        label_map[old_label] = i


    new_mods = []
    for label in mods:
        if label == -1:
            new_mods.append(-1)   
        else:
            new_mods.append(label_map[label])
    mods = np.array(new_mods, dtype=int)

    return mods

def analyze_module(mods, module_num, data_main, adSC, out_dir):
    module_indices = np.where(mods == module_num)[0]
    A = adSC['Adj'][np.ix_(module_indices, module_indices)]
    G = nx.from_numpy_array(A) 
    G_Distance  = nx.from_numpy_array(1 - A)  

    deg = dict(G.degree(weight='weight'))     
    clo = nx.closeness_centrality(G_Distance, distance='weight') 
    bet = nx.betweenness_centrality(G_Distance, weight='weight') 
    eig = nx.eigenvector_centrality_numpy(G, weight='weight') 

    df = pd.DataFrame({
        'gene_names': [data_main['name'][i] for i in module_indices],
        'ext_names': [data_main['extName'][i] for i in module_indices],
        'degree': pd.Series(deg),
        'closeness': pd.Series(clo),
        'betweenness': pd.Series(bet),
        'eigenvector': pd.Series(eig)
    })

    for col in ['degree', 'closeness', 'betweenness', 'eigenvector']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df['degRank']   = df['degree'].rank(ascending=False, method='dense')
    df['closeRank'] = df['closeness'].rank(ascending=False, method='dense')
    df['eigRank']   = df['eigenvector'].rank(ascending=False, method='dense')
    df['betRank']   = df['betweenness'].rank(ascending=False, method='dense')
    df['avgRank']   = df[['degRank', 'closeRank', 'eigRank','betRank']].mean(axis=1)

    df.to_csv(f"{out_dir}/module_{module_num}.csv", index=False, float_format='%.4f')
    return df


def write_compare_modules(csv_directory):
    files = glob.glob(f'{csv_directory}/module_*.csv')
    summary = {
        'identifier': [], 'size': [], 'medianDegree': [],
        'maxDegree': [], 'correlatedCondition': [], 'hubGene': []
    }

    for f in files:
        df = pd.read_csv(f)
        mod_num = int(f.split('_')[-1].split('.')[0])
        df['avgRank'] = df[['degRank', 'closeRank', 'eigRank','betRank']].mean(axis=1)
        hub = df.loc[df['avgRank'].idxmin(), 'ext_names']
        summary['identifier'].append(mod_num)
        summary['size'].append(len(df))
        summary['medianDegree'].append(df['degree'].median())
        summary['maxDegree'].append(df['degree'].max())
        summary['correlatedCondition'].append('N/A')
        summary['hubGene'].append(hub)

    return pd.DataFrame(summary)


def plot_to_bytes(fig) -> bytes:
    """
    Serialize a Matplotlib figure to bytes for download.
    """
    buf = BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    return buf.getvalue()
