#Importing libraries
import matplotlib.pyplot as plt
import seaborn as sns
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import pdist, squareform
import umap
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from itertools import product
import warnings
import networkx as nx

#Importing the Files
rlogCounts = pd.read_csv('deseq2_rlog_normalized_counts.txt', delimiter = "\t", index_col = 0)
counts = pd.read_csv("deseq2_raw_counts.txt", delimiter="\t", index_col=0)
sampleInfo = pd.read_csv("samplesTableSeedless.csv")

def columnFilter_loc(df_to_filter, ref_df, targ_col, loc_cond, target_col_toExport):
    if isinstance(target_col_toExport, str):
        target_col_toExport = [target_col_toExport]
    filtered_df = ref_df.loc[ref_df[targ_col] == loc_cond, target_col_toExport]
    matching_col = df_to_filter.columns.isin(filtered_df.values.ravel())
    return df_to_filter.loc[:, matching_col]



name = counts.index.astype(str).tolist() 
nc1 = columnFilter_loc(rlogCounts, sampleInfo, 'source', 'far', 'sample')  # Normalized Count for 'far'
nc2 = columnFilter_loc(rlogCounts, sampleInfo, 'source', 'near', 'sample')  # Normalized Count for 'near'
rc1 = columnFilter_loc(counts, sampleInfo, 'source', 'far', 'sample')  # Raw Count for 'far'
rc2 = columnFilter_loc(counts, sampleInfo, 'source', 'near', 'sample')  # Raw Count for 'near'


data = {
    "name": name,
    "extName": name,
    "nc1": nc1,
    "rc1": rc1,
    "nc2": nc2,
    "rc2": rc2
}



print('length of extName is:', len(data['extName']))
print('Class of extName is:', type(data['extName']))
print(data['extName'])


batch_size = 1000
url = "https://rest.ensembl.org/lookup/id"
headers = {"Content-Type": "application/json"}


def fetch_gene_data_batch(batch):
    data = {"ids": batch}
    response = requests.post(url, json=data, headers=headers)

    results = []
    if response.status_code == 200:
        gene_infos = response.json()
        for gene_id in batch:
            gene_info = gene_infos.get(gene_id)
            if gene_info:
                length = abs(gene_info.get("end", 0) - gene_info.get("start", 0))
                results.append({
                    "id": gene_id,
                    "hasEnsData": True,
                    "name": gene_info.get("display_name", ""),
                    "type": gene_info.get("biotype", ""),
                    "length": length
                })
            else:
                results.append({
                    "id": gene_id,
                    "hasEnsData": False,
                    "name": "",
                    "type": "",
                    "length": ""
                })
    else:
        for gene_id in batch:
            results.append({
                "id": gene_id,
                "hasEnsData": False,
                "name": "",
                "type": "",
                "length": ""
            })

    return results

def split_batches(gene_ids, batch_size):
    for i in range(0, len(gene_ids), batch_size):
        yield gene_ids[i:i + batch_size]

def fetch_data_concurrently(gene_ids, batch_size, max_workers=10):
    futures = []
    with ThreadPoolExecutor(max_workers = max_workers) as executor:
        for batch in split_batches(gene_ids, batch_size):
            future = executor.submit(fetch_gene_data_batch, batch)
            futures.append((batch, future))

    results = []
    for batch, future in futures:
        results.extend(future.result())

    return results

gene_data = fetch_data_concurrently(name, batch_size)

gene_df = pd.DataFrame(gene_data)

print(gene_df.head())
gene_df.to_csv("ensembl_gene_info_batch_filtered.csv", index=False)


def update_extName(data, gene_df):
    gene_map = dict(zip(gene_df['id'], gene_df['name']))
    for i, gene_id in enumerate(data['extName']):
        if gene_id in gene_map:
            data['extName'][i] = gene_map[gene_id]
    return data

updated_data = update_extName(data, gene_df)

print(len(updated_data['name']))


updated_data_new = {
    "name": gene_df['id'].tolist(),
    "extName": updated_data['extName'],
    "nc1": nc1,
    "rc1": rc1,
    "nc2": nc2,
    "rc2": rc2
}


#print(updated_data_new['name'])
print(len(updated_data_new['name']))
#print(updated_data_new['extName'])
print(len(updated_data_new['extName']))

print(len(updated_data_new['rc2']))
print(len(updated_data_new['rc1']))
print(len(updated_data_new['nc1']))
print(len(updated_data_new['nc2']))

min_rc = 100
min_nc = 0
updated_data_new['rc1'] = updated_data_new['rc1'].apply(pd.to_numeric, errors='coerce')
updated_data_new['rc2'] = updated_data_new['rc2'].apply(pd.to_numeric, errors='coerce')
updated_data_new['nc1'] = updated_data_new['nc1'].apply(pd.to_numeric, errors='coerce')
updated_data_new['nc2'] = updated_data_new['nc2'].apply(pd.to_numeric, errors='coerce')
total_raw_counts = updated_data_new['rc1'].sum(axis=1, skipna=True) + updated_data_new['rc2'].sum(axis=1, skipna=True)
total_normalized_counts = updated_data_new['nc1'].sum(axis=1, skipna=True) + updated_data_new['nc2'].sum(axis=1, skipna=True)
high_raw_counts = (total_raw_counts >= min_rc) & (total_normalized_counts > min_nc)
print(high_raw_counts)
print(sum(high_raw_counts))



print(total_raw_counts)
print(total_normalized_counts)
less_than_0_normalized = total_normalized_counts[total_normalized_counts < 0].count()
less_than_100_raw = total_raw_counts[total_raw_counts < 100].count()

greater_than_0_normalized = total_normalized_counts[total_normalized_counts >= 0].count()
greater_than_100_raw = total_raw_counts[total_raw_counts >= 100].count()

print(f"Number of raw values less than 100: {less_than_100_raw}")
print(f"Number of norm values less than 0: {less_than_0_normalized}")
print(f"Number of raw values Greater than or equal to 100: {greater_than_100_raw}")
print(f"Number of Norm values Greater than or equal to 0: {greater_than_0_normalized}")
print(f"Number of total gene for Normalized: {greater_than_0_normalized+ less_than_0_normalized}")
print(f"Number of total gene for Normalized: {greater_than_100_raw+ less_than_100_raw}")


print(f"{((total_raw_counts >= 100) & (total_normalized_counts >= 0)).sum()}")
print(f"{(17317 - ((total_raw_counts < 100) & (total_normalized_counts < 0)).sum())}")
print(f"{((total_raw_counts < 100) & (total_normalized_counts >= 0)).sum()}")
print(f"{(17317 - ((total_raw_counts < 100) & (total_normalized_counts > 0)).sum())}")
print(f"{(((total_raw_counts >= 100) & (total_normalized_counts < 0)).sum())}")







def filter_data_by_iteration(data, high_raw_counts):


    filtered_name = []
    filtered_extName = []
    filtered_rc1 = []
    filtered_rc2 = []
    filtered_nc1 = []
    filtered_nc2 = []


    for i, gene_id in enumerate(data['name']):
        if high_raw_counts.get(gene_id, True):  # If True, keep the data
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


filtered_data_1count = filter_data_by_iteration(updated_data_new, high_raw_counts)


min_samp = 3
s1 = (updated_data_new['rc1'] > 0).sum(axis=1)
s2 = (updated_data_new['rc2'] > 0).sum(axis=1)

more_samples = (s1 >= min_samp) & (s2 >= min_samp)
print(more_samples)


s1_jittered = s1 + np.random.uniform(-0.2, 0.2, size=s1.shape)
s2_jittered = s2 + np.random.uniform(-0.2, 0.2, size=s2.shape)


plt.scatter(s1_jittered, s2_jittered, s=10, alpha=0.5)


plt.xlabel('Far Samples')
plt.ylabel('Near Samples')

plt.xlim(-0.5, 18.5)
plt.ylim(-0.5, 18.5)


plt.xticks(range(0, 19))
plt.yticks(range(0, 19))
plt.title('Scatter Plot of Sample Counts')


plt.show()



min_samp = 3

def generate_non_zero_mask(data, ids, min_samp):


    rc1_df = pd.DataFrame(data['rc1'])
    rc2_df = pd.DataFrame(data['rc2'])

    s1 = (rc1_df > 0).sum(axis=1)
    s2 = (rc2_df > 0).sum(axis=1)

    more_samples = (s1 >= min_samp) & (s2 >= min_samp)

    mask_series = pd.Series(more_samples.values, index=ids)
    mask_series.index.name = 'id'
    return mask_series


fewsample_mask = generate_non_zero_mask(filtered_data_1count, filtered_data_1count['name'], min_samp)
print(fewsample_mask)


empty_extName_count = sum(1 for extName in filtered_data_1count['extName'] if extName.strip() == '')

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


filtered_data_2moresample = filter_data_by_moresample(filtered_data_1count, fewsample_mask)
print(len(filtered_data_2moresample['name']))

def filter_extName_empty_or_space(data):


    filtered_name = []
    filtered_extName = []
    filtered_rc1 = []
    filtered_rc2 = []
    filtered_nc1 = []
    filtered_nc2 = []


    for i, ext_name in enumerate(data['extName']):
        if ext_name.strip():
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

filtered_data_3noname = filter_extName_empty_or_space(filtered_data_2moresample)


def filter_and_compare_gene_df(data, gene_df, min_len):



    gene_df['length'] = pd.to_numeric(gene_df['length'], errors='coerce')


    not_coding_gene = gene_df['type'] != 'protein_coding'
    too_short = gene_df['length'] < min_len
    no_data = ~gene_df['hasEnsData']


    mask = no_data | not_coding_gene | too_short


    filtered_gene_names = gene_df.loc[mask, 'name']


    filtered_extName = [ext_name for ext_name in data['extName'] if ext_name not in filtered_gene_names.tolist()]


    filtered_indices = [i for i, ext_name in enumerate(data['extName']) if ext_name not in filtered_gene_names.tolist()]

    filtered_data = {
        'name': [data['name'][i] for i in filtered_indices],
        'extName': [data['extName'][i] for i in filtered_indices],
        'rc1': [data['rc1'][i] for i in filtered_indices],
        'rc2': [data['rc2'][i] for i in filtered_indices],
        'nc1': [data['nc1'][i] for i in filtered_indices],
        'nc2': [data['nc2'][i] for i in filtered_indices]
    }

    return filtered_data, filtered_gene_names


min_len = 100
filtered_data_4_otherfilters, filtered_gene_names_final = filter_and_compare_gene_df(filtered_data_3noname, gene_df, min_len)


print("Filtered gene names_final:", filtered_gene_names_final)

print(filtered_data_4_otherfilters)


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

results = compute_correlations(filtered_data_4_otherfilters)

# Print the results
print("FDR 10% P-value:", results['FDR_thresholds'][0.1]['pstar'])
print("FDR 10% Correlation:", results['FDR_thresholds'][0.1]['cstar'])
print("FDR 50% P-value:", results['FDR_thresholds'][0.5]['pstar'])
print("FDR 50% Correlation:", results['FDR_thresholds'][0.5]['cstar'])
print("FDR 90% P-value:", results['FDR_thresholds'][0.9]['pstar'])
print("FDR 90% Correlation:", results['FDR_thresholds'][0.9]['cstar'])



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def find_strongly_connected_genes(corrs, cstar10, cstar50, min_corr_threshold=200):
    

    sum1 = np.sum(corrs['C1'] > cstar10, axis=1) + np.sum(corrs['C2'] > cstar10, axis=1)
    sum2 = np.sum(corrs['C1'] > cstar50, axis=1) + np.sum(corrs['C2'] > cstar50, axis=1)

    plt.hist(sum1, bins=100)
    plt.title("Number of correlations above 10% FDR")
    plt.xlabel("Count")
    plt.ylabel("Frequency")
    plt.show()

    plt.hist(sum2, bins=100)
    plt.title("Number of correlations above 50% FDR")
    plt.xlabel("Count")
    plt.ylabel("Frequency")
    plt.show()

    plt.scatter(sum2, sum1, s=1)
    plt.xlabel("Number of correlations below 50% FDR")
    plt.ylabel("Number of correlations below 10% FDR")
    plt.title("Number of strong correlations (each point is a gene)")
    plt.show()

    testSC = sum2 > min_corr_threshold
    return testSC

def find_strongly_connected_genes_subset(corrs, test):
    
    corrsSC = {}
    corrsSC['C1'] = corrs['C1'][np.ix_(test, test)]
    corrsSC['C2'] = corrs['C2'][np.ix_(test, test)]
    return corrsSC


corrs = {
    'C1': results['C1_corr'], 
    'C2': results['C2_corr']   
}

cstar10 = results['FDR_thresholds'][0.1]['cstar']
cstar50 = results['FDR_thresholds'][0.5]['cstar']


testSC = find_strongly_connected_genes(corrs, cstar10, cstar50)
corrsSC = find_strongly_connected_genes_subset(corrs, testSC)

print("Number of strongly connected genes:", np.sum(testSC))



def compute_adjacency(corrs, cap=1):

    x1 = np.copy(corrs['C1'])
    x2 = np.copy(corrs['C2'])

    x1[x1 < 0.001] = 0
    x2[x2 < 0.001] = 0


    A = np.abs(x1 - x2) 
    A[A > cap] = cap  
    A /= cap  
    np.fill_diagonal(A, 1)  

    plt.hist(np.log(A[A > 0]).flatten(), bins=100)
    plt.title("Histogram of log-transformed adjacency values")
    plt.xlabel("Log(adj)")
    plt.ylabel("Frequency")
    plt.show()


    T = pairwise_distances(A, metric='euclidean')

    return {
        'cap': cap,
        'Adj': A,
        'Tom': T
    }

adSC = compute_adjacency(corrsSC, cap=1)



import umap
import pandas as pd

def quick_umap(T, n_neighbors=15, min_dist=0.1, seed=123, n_components=2):

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        random_state=seed,
        metric='euclidean'
    )
    umap_embedding = reducer.fit_transform(T)
    if n_components == 2:
        return pd.DataFrame({'x': umap_embedding[:, 0], 'y': umap_embedding[:, 1]})
    else:
        return pd.DataFrame({'x': umap_embedding[:, 0], 'y': umap_embedding[:, 1], 'z': umap_embedding[:, 2]})

# Create UMAP embedding
umapSC = quick_umap(adSC['Tom'], n_neighbors=30, min_dist=0.1, seed=123, n_components=2)

# Plot the UMAP embedding
plt.scatter(umapSC['x'], umapSC['y'], s=4, color='blue')
plt.title('UMAP of Strongly Connected Genes')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.show()




def density_cluster(umap_embedding, eps, min_samples, min_mod_size, merge0=False):
    
    dbscan = DBSCAN(eps=eps, min_samples=min_samples).fit(umap_embedding)
    mods = dbscan.labels_


    for m in range(max(mods) + 1):
        if np.sum(mods == m) <= min_mod_size:
            mods[mods == m] = -1

   
    if merge0:
        noise_indices = np.where(mods == -1)[0]
        for idx in noise_indices:

            distances = np.linalg.norm(umap_embedding.values - umap_embedding.iloc[idx].values, axis=1)
            nearest_cluster_idx = np.argmin(distances)
            mods[idx] = mods[nearest_cluster_idx]

    return mods


#modsSC = density_cluster(umapSC, eps=0.55, min_samples=18, min_mod_size=50, merge0=True)
modsSC = density_cluster(umapSC, eps=0.7, min_samples=20, min_mod_size=50, merge0=False)

plt.scatter(umapSC['x'], umapSC['y'], c=modsSC, cmap='viridis', s=4)
plt.title('Density-Based Clustering of Strongly Connected Genes')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.colorbar(label='Cluster Label')
plt.show()




def evaluate_clusters(data, labels):
    
    metrics = {}
    if len(set(labels)) > 1:  
        metrics['Silhouette Score'] = silhouette_score(data, labels)
        metrics['Davies-Bouldin Index'] = davies_bouldin_score(data, labels)
        metrics['Calinski-Harabasz Index'] = calinski_harabasz_score(data, labels)
    else:
        metrics['Silhouette Score'] = None
        metrics['Davies-Bouldin Index'] = None
        metrics['Calinski-Harabasz Index'] = None

    return metrics


cluster_metrics = evaluate_clusters(umapSC, modsSC)
print(cluster_metrics)


#This is the working code
import umap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score, davies_bouldin_score
from itertools import product
import warnings


warnings.filterwarnings('ignore')


def quick_umap(T, n_neighbors, min_dist, n_components, metric, seed):
    
    T = T.astype(np.float32)  
    reducer = umap.UMAP(
        n_neighbors=int(n_neighbors),
        min_dist=float(min_dist),
        n_components=n_components,
        metric=metric,
        random_state=seed
    )
    return reducer.fit_transform(T)

def density_cluster(umap_embedding, eps, min_samples):
    
    dbscan = DBSCAN(eps=eps, min_samples=min_samples).fit(umap_embedding)
    return dbscan.labels_


def grid_search_umap_dbscan(data, n_neighbors_range, min_dist_range, eps_range, min_samples_range, metric='euclidean', seed=123):
    
    results = []
    for n_neighbors, min_dist, eps, min_samples in product(n_neighbors_range, min_dist_range, eps_range, min_samples_range):
        try:
            
            umap_embedding = quick_umap(data, n_neighbors, min_dist, n_components=2, metric=metric, seed=seed)

            
            labels = density_cluster(umap_embedding, eps, min_samples)

            
            if len(set(labels)) > 1:
                silhouette_avg = silhouette_score(umap_embedding, labels)
                davies_bouldin = davies_bouldin_score(umap_embedding, labels)
            else:
                silhouette_avg = -1
                davies_bouldin = np.inf

            # Store results
            results.append({
                'n_neighbors': n_neighbors,
                'min_dist': min_dist,
                'eps': eps,
                'min_samples': int(min_samples),
                'silhouette_score': silhouette_avg,
                'davies_bouldin': davies_bouldin
            })
        except Exception as e:
            print(f"Error at n_neighbors={n_neighbors}, min_dist={min_dist}, eps={eps}, min_samples={min_samples}: {e}")
    return pd.DataFrame(results)


if __name__ == "__main__":
    
    np.random.seed(123)
    data =  adSC['Tom']

    
    data = np.nan_to_num(data).astype(np.float32)

    
    n_neighbors_range = [10, 20, 30]
    min_dist_range = [0.1, 0.2, 0.3]
    eps_range = [0.3, 0.5, 0.7]
    min_samples_range = [5, 10, 15]

    
    results = grid_search_umap_dbscan(data, n_neighbors_range, min_dist_range, eps_range, min_samples_range)

    
    results = results.sort_values(by='silhouette_score', ascending=False)
    print("Top Grid Search Results:")
    print(results.head())

    
    best_params = results.iloc[0]
    print(f"Best Parameters: {best_params}")

    
    optimal_umap = quick_umap(
        data,
        n_neighbors=best_params['n_neighbors'],
        min_dist=best_params['min_dist'],
        n_components=2,
        metric='euclidean',
        seed=123
    )
    optimal_labels = density_cluster(optimal_umap, eps=best_params['eps'], min_samples=int(best_params['min_samples']))
    
    optimal_labels = optimal_labels.astype(int)


    cmap = plt.cm.get_cmap('viridis', len(np.unique(optimal_labels)))


results_df = pd.DataFrame(results)

# Save results to a CSV file
results_df.to_csv("umap_dbscan_grid_search_results.csv", index=False)
print("Grid search results saved to 'umap_dbscan_grid_search_results.csv'")

#Finding out the List of top Parameters


results = pd.read_csv('umap_dbscan_grid_search_results.csv')
# Normalize Silhouette Score (maximize)
results['silhouette_norm'] = (results['silhouette_score'] - results['silhouette_score'].min()) / (
    results['silhouette_score'].max() - results['silhouette_score'].min()
)

# Normalize DBI (minimize)
results['dbi_norm'] = 1 - (results['davies_bouldin_score'] - results['davies_bouldin_score'].min()) / (
    results['davies_bouldin_score'].max() - results['davies_bouldin_score'].min()
)

# Composite Score
results['composite_score'] = 0.3*results['silhouette_norm'] + 0.7*results['dbi_norm']

# Sort by Composite Score
sorted_results = results.sort_values(by='composite_score', ascending=False)

# Display the top 10 parameters
top_results = sorted_results.head(15)
top_results.to_csv("top_clustering_parameters.csv", index=False)
print("Top clustering parameters saved to 'top_clustering_parameters.csv'")


#Plot based on top clustering parameters
def quick_umap(T, n_neighbors, min_dist, n_components, metric, seed):
   
    T = T.astype(np.float32)  
    reducer = umap.UMAP(
        n_neighbors=int(n_neighbors),
        min_dist=float(min_dist),
        n_components=n_components,
        metric=metric,
        random_state=seed
    )
    return reducer.fit_transform(T)

def plot_umap_dbscan(umap_embedding, labels, title, silhouette, davies_bouldin, ax):
  
    labels = labels.astype(int)
    cmap = plt.cm.get_cmap('viridis', len(np.unique(labels)))  
    scatter = ax.scatter(umap_embedding[:, 0], umap_embedding[:, 1], c=labels, cmap=cmap, s=10)
    ax.set_title(f"{title}\nSilhouette: {silhouette:.2f}, DBI: {davies_bouldin:.2f}")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")


top_results_file = "top_clustering_parameters.csv"  
top_results = pd.read_csv(top_results_file)


top_results = top_results.sort_values(by='composite_score', ascending=False)


num_plots = min(12, len(top_results)) 
fig, axes = plt.subplots(4, 3, figsize=(18, 16)) 
axes = axes.flatten()


for i, (_, res) in enumerate(top_results.iterrows()):
    if i >= num_plots:
        break

    umap_embedding = quick_umap(
        data, 
        n_neighbors=int(res['n_neighbors']), 
        min_dist=res['min_dist'], 
        n_components=2, 
        #metric=res['metric'],
        metric = 'euclidean',
        seed=123 + i  
    )

 
    labels = density_cluster(umap_embedding, eps=res['eps'], min_samples=int(res['min_samples']))

  
    plot_umap_dbscan(
        umap_embedding,
        labels,
        title=f"UMAP(n_neighbors={res['n_neighbors']}, min_dist={res['min_dist']}, metric={res['metric']})\n"
              f"DBSCAN(eps={res['eps']}, min_samples={int(res['min_samples'])})",
        
        silhouette=res['silhouette_score'],
        davies_bouldin=res['davies_bouldin_score'],
        ax=axes[i]
    )


for j in range(num_plots, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.show()



def analyze_module(mods, module_num, data, adSC, corrsSC):

    module_indices = np.where(mods == module_num)[0]
    adj_matrix = adSC['Adj'][np.ix_(module_indices, module_indices)]


    G = nx.from_numpy_array(adj_matrix)


    degrees = dict(G.degree(weight='weight'))
    closeness = nx.closeness_centrality(G)
    betweenness = nx.betweenness_centrality(G, weight='weight')
    eigenvector = nx.eigenvector_centrality_numpy(G, weight='weight')

    module_data = {
        'size': len(module_indices),
        'degree': list(degrees.values()),
        'closeness': list(closeness.values()),
        'betweenness': list(betweenness.values()),
        'eigenvector': list(eigenvector.values()),
        'gene_names': [data['name'][i] for i in module_indices],
        'ext_names': [data['extName'][i] for i in module_indices]
    }


    df = pd.DataFrame(module_data)
    df.to_csv(f'module_{module_num}.csv', index=False)

    return module_data


mod_list = []
for m in np.unique(modsSC):
    if m != -1:
        mod_list.append(analyze_module(modsSC, m, filtered_data_4_otherfilters, adSC, corrsSC))


import pandas as pd
import numpy as np
import glob

def write_compare_modules(csv_directory, output_file):


    module_files = glob.glob(f'{csv_directory}/module_*.csv')

    summary_data = {
        'identifier': [],
        'size': [],
        'medianDegree': [],
        'maxDegree': [],
        'correlatedCondition': [],
        'hubGene': []
    }

    for file in module_files:

        module_num = int(file.split('_')[-1].split('.')[0])


        df = pd.read_csv(file)


        size = df.shape[0]


        median_degree = df['degree'].median()
        max_degree = df['degree'].max()


        avg_rank = df[['degree', 'closeness', 'eigenvector']].rank(ascending=False).mean(axis=1)


        hub_gene_idx = avg_rank.argmin()
        hub_gene_name = df.loc[hub_gene_idx, 'gene_names']


        summary_data['identifier'].append(module_num)
        summary_data['size'].append(size)
        summary_data['medianDegree'].append(median_degree)
        summary_data['maxDegree'].append(max_degree)
        summary_data['correlatedCondition'].append('N/A')
        summary_data['hubGene'].append(hub_gene_name)


    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_file, index=False)


write_compare_modules(csv_directory='.', output_file='./compareModules_nocap_2.14.csv')



import pandas as pd


module_df = pd.read_csv('module_2.csv')


numeric_cols = ['degree', 'closeness', 'eigenvector', 'betweenness']
module_df[numeric_cols] = module_df[numeric_cols].apply(pd.to_numeric, errors='coerce')


module_df['degRank'] = module_df['degree'].rank(ascending=False, method='dense').astype(float)
module_df['closeRank'] = module_df['closeness'].rank(ascending=False, method='dense').astype(float)
module_df['eigRank'] = module_df['eigenvector'].rank(ascending=False, method='dense').astype(float)
module_df['betRank'] = module_df['betweenness'].rank(ascending=False, method='dense').astype(float)


module_df['avgRank'] = module_df[['degRank', 'closeRank', 'eigRank']].mean(axis=1)


rank_cols = ['degRank', 'closeRank', 'eigRank', 'avgRank', 'betRank']
module_df[rank_cols] = module_df[rank_cols].apply(pd.to_numeric)


module_df.to_csv('module_2_updated_1.21.25.csv', index=False, float_format='%.2f')

print("The ranks have been added successfully and saved as 'module_2_updated.csv'.")