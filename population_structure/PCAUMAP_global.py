import umap.umap_ as umap
# from scipy.sparse.csgraph import connected_components
import sklearn.cluster as cluster
import hdbscan, warnings
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd

df = pd.read_csv("genotype_HpGP2_decomposed_prunned.txt", sep="\t", index_col=0)
df = df.T

warnings.resetwarnings()

np.random.seed(seed=42)

n_neighbors = 100
n_components = 2
n_clusters = 5

model_pca = PCA(random_state=42, n_components=n_components, svd_solver='arpack')
pca_x = model_pca.fit_transform(df)

clusterable_embedding = umap.UMAP(init='random', min_dist=0, n_neighbors = n_neighbors, n_components=n_components, random_state=42)

vecs_list = clusterable_embedding.fit_transform(pca_x)


kmeans_labels = sklearn.cluster.KMeans(init='k-means++', n_clusters=n_clusters, random_state=42).fit_predict(vecs_list)
spectral_labels = kmeans_labels+1

fig = plt.figure(figsize=(20.0, 10.0))

x = vecs_list[:,0]
y = vecs_list[:,1]

cm = plt.cm.get_cmap('tab20')

marker_list = ['o', 'D', 'v', 's', '*', '+' ]

cluster_names = ["Latin America (LA)", "East Asia (EA)", "Africa (AF)", "Europe (EU)", "Asia-Amerind (AA)"]

for i in range(1, n_clusters+1):
    x_draw = x[spectral_labels == i]
    y_draw = y[spectral_labels == i]
    mod = (i-1)%len(marker_list)
    cluster_name = "Global"+str(i)
    plt.scatter(x_draw, y_draw, s=50, label=cluster_names[i-1], cmap='jet', marker=marker_list[mod], alpha=0.7)
    

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=18)
plt.tight_layout()
name = "global_"+str(n_components)+"_"+str(n_neighbors)+"_"+str(n_clusters)

fig.savefig(name+".png")
fig.show()
