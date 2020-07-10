import velocyto as vcy
import glob
from singlecellexperiment import SingleCellExperiment
import loompy
from sklearn.neighbors import NearestNeighbors
import igraph
import numpy as np
import matplotlib.pyplot as plt
import sys

rdata = eval(sys.argv[1])
celltype = sys.argv[2]
output = sys.argv[3]
svg = sys.argv[4]
looms = eval(sys.argv[5])

reduction = "UMAPCELLTYPE"

print("Loading SCE...")
sce = SingleCellExperiment.fromRData(rdata[celltype])
looms = eval(looms)
looms = glob.glob("/juno/work/shah/ceglian/rnascp/velocity/velocity*.loom")
loompy.combine(list(looms), output, key="Accession")
print("Loading Combined Loom...")
vlm = vcy.VelocytoLoom(output)

print("Filtering to SCE...")
unfiltered_embedding = sce.getReducedDims(reduction)
barcodes = sce.colnames
valid_barcodes = [barcode.split("_")[1] for barcode in barcodes] 
unfiltered_barcodes = [x.split(":")[1].rstrip("x") for x in vlm.ca["CellID"]]
visited = set()
filtered_barcodes = []
for barcode in unfiltered_barcodes:
    if barcode in valid_barcodes and barcode not in visited:
        filtered_barcodes.append(True)
        visited.add(barcode)
    else:
        filtered_barcodes.append(False)

print("Filtering Cells...")
vlm.filter_cells(filtered_barcodes)

print("Scoring genes...")
vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
vlm.score_cv_vs_mean(3000, plot=False, max_expr_avg=35)

print("Normalization...")
vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())

print("Setting Clusters and Filtering Embedding...")
clusters = []
embedding = []
unfiltered_barcodes = [x.split(":")[1].rstrip("x") for x in vlm.ca["CellID"]]
cluster_ids = sce.colData["seurat_clusters"]
visited = set()
for i in range(len(cluster_ids)):
    if valid_barcodes[i] in unfiltered_barcodes and valid_barcodes[i] not in visited:
        clusters.append(cluster_ids[i])
        embedding.append(unfiltered_embedding[i])
        visited.add(valid_barcodes[i])

assert len(clusters) == len(vlm.ca["CellID"]), "{}-{}".format(len(clusters), len(vlm.ca["CellID"])) 
vlm.set_clusters(clusters)

print("Scoring Clusters...")
vlm.score_cluster_expression(min_avg_U=0.007, min_avg_S=0.06)
print("Filtering genes...")
vlm.filter_genes(by_detection_levels=True, by_cluster_expression=True,by_cv_vs_mean=True)

print("Adjusting Splice/Unspliced...")
vlm.normalize_by_total(min_perc_U=0.5)
vlm.adjust_totS_totU(normalize_total=True, fit_with_low_U=False, svr_C=1, svr_gamma=1e-04)

print("Running PCA...")
vlm.perform_PCA()
n_comps = np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_))>0.0055))[0][0]
vlm.pcs[:,1] *= -1

print("Creating KNN Graph...")
nn = NearestNeighbors(50)
nn.fit(vlm.pcs[:,:4])
knn_pca = nn.kneighbors_graph(mode='distance')
knn_pca = knn_pca.tocoo()
G = igraph.Graph(list(zip(knn_pca.row, knn_pca.col)), directed=False, edge_attrs={'weight': knn_pca.data})
VxCl = G.community_multilevel(return_levels=False, weights="weight")
labels = np.array(VxCl.membership)

print("Annotating Clusters...")
manual_annotation = {str(i):[i] for i in labels}
annotation_dict = {v:k for k, values in manual_annotation.items() for v in values }
clusters = np.array([annotation_dict[i] for i in labels])
colors20 = np.vstack((plt.cm.tab20b(np.linspace(0., 1, 20))[::2], plt.cm.tab20c(np.linspace(0, 1, 20))[1::2]))
vlm.set_clusters(clusters, cluster_colors_dict={k:colors20[v[0] % 20,:] for k,v in manual_annotation.items()})

print("Imputing...")
vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=16)

print("Fitting Gamma...")
vlm.normalize_median()
vlm.fit_gammas()

print("Normalizing Imputation...")
vlm.normalize(which="imputed", size=False, log=True)
vlm.Pcs = np.array(vlm.pcs[:,:2], order="C")

print("Calculating Velocity...")
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift()

vlm.ts = sce.getReducedDims("umapcelltype")
print("Extrapolating Transitions...")
vlm.extrapolate_cell_at_t(delta_t=1)
vlm.estimate_transition_prob(hidim="Sx_sz", embed="Pcs", transform="log", psc=1,
                             n_neighbors=150, knn_random=True, sampled_fraction=1)

vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

print("Calculating Vectors...")
vlm.calculate_grid_arrows(smooth=0.9, steps=(25, 25), n_neighbors=200)

print("Plotting...")
plt.figure(None,(9,9))
vlm.plot_grid_arrows(scatter_kwargs_dict={"alpha":0.7, "lw":0.7, "edgecolor":"0.4", "s":70, "rasterized":True},
                     min_mass=2.9, angles='xy', scale_units='xy',
                     headaxislength=2.75, headlength=5, headwidth=4.8, quiver_scale=0.35, scale_type="absolute")

plt.gca().invert_xaxis()
plt.axis("off")
plt.axis("equal");
plt.title("{} Trajectory".format(celltype))
plt.savefig(svg)
print("Complete!")

