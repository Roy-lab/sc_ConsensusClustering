"""
 Written by Junha Shin (junha.shin@wisc.edu)
 Copyright (c) 2021 by Junha Shin and Sushmita Roy, 
 Wisconsin Institute for Discovery, University of Wisconsin-Madison
 All rights reserved.
"""
import sys, re, os
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial import distance_matrix
import igraph as ig
import louvain

matfile = sys.argv[1]	# data matrix file, rows will be clustered
k = int(sys.argv[2])	# knn k neighborhood number
cl = int(sys.argv[3])	# clustering k number
outname = sys.argv[4]	# output naming
# Results
rmat = pd.DataFrame()
resT = pd.DataFrame(columns=['resolution'])
rmat_name = 'louvain_cluster_' + outname + '.txt'
resT_name = 'used_resolution_' + outname + '.txt'

def knn_bin_adjmat(smat,k):
	knnsim = np.zeros((len(smat),len(smat)))
	row = 0
	for r in smat:
		ids = np.flip(np.argsort(r))    # descending sort
		for col in ids[0:k+1]:
			knnsim[row,col] = 1
			knnsim[col,row] = 1
		row = row+1
	return knnsim

def get_igraph_from_adjacency(adjacency, directed=None):
	sources, targets = adjacency.nonzero()
	weights = adjacency[sources, targets]
	if isinstance(weights, np.matrix):
		weights = weights.A1
	g = ig.Graph(directed=directed)
	g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
	g.add_edges(list(zip(sources, targets)))
	try:
		g.es['weight'] = weights
	except:
		pass
	return g

def getNClusters(graph,n_cluster,range_min=0,range_max=3,max_steps=20):
# adopted from: https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018_bulkpeaks/run_clustering_buenrostro2018bulkpeaks.ipynb
	this_step = 0
	this_min = float(range_min)
	this_max = float(range_max)
	while this_step < max_steps:
		print('step ' + str(this_step))
		this_resolution = this_min + ((this_max-this_min)/2)
		# louvain clustering part
		partition_type = louvain.RBConfigurationVertexPartition
		partition_kwargs = dict()
		partition_kwargs["resolution_parameter"] = this_resolution
		weights = np.array(G.es["weight"]).astype(np.float64)
		partition_kwargs["weights"] = weights
		part = louvain.find_partition(G,partition_type,**partition_kwargs)
		clusters = np.array(part.membership)
		this_clusters = len(np.unique(clusters))
		print('got ' + str(this_clusters) + ' at resolution ' + str(this_resolution))
		if this_clusters > n_cluster:
			this_max = this_resolution
		elif this_clusters < n_cluster:
			this_min = this_resolution
		else:
			return(clusters+1, this_resolution)     # +1 for clusterIDs
		this_step += 1
	print('Cannot find the number of clusters')
	print('Clustering solution from last iteration is used:' + str(this_clusters) + ' at resolution ' + str(this_resolution))
	return [], 0

mat = pd.read_csv(matfile,sep='\t',header=0,index_col=0)
print ('Input matrix shape=' + str(mat.values.shape))

# Euclidean distance matrix
print ('Calculating distance matrix')
eudistmat = distance_matrix(mat.values, mat.values)
# similarity matrix normalization by Gaussian kernel
sig = np.std(eudistmat)
print ('distance matrix stdev=' + str(sig))
simmat = np.exp(-((eudistmat)/(2*sig*sig)))

print ('kNN=' + str(k))
knn_binmat = knn_bin_adjmat(simmat,k)
simdistmat = 1-simmat
knnmat=np.multiply(simdistmat,knn_binmat)		# similarity=1 is assigned as "0" in this matrix
simdistmin=np.min(simdistmat[np.nonzero(simdistmat)])   # non-zero min distance
simdist_smallval = simdistmin * 0.1     		# multiplying 0.1 to assign most small distance
print ('smallest distance=' + str(simdist_smallval))
smallvalidx=np.multiply((simdistmat==0),(knn_binmat==1))# index of similarity=1 edges
smallvalmat=np.multiply(simdist_smallval,smallvalidx)   # matrix of similarity=1 edges
knnmat = knnmat + smallvalmat

G = get_igraph_from_adjacency(np.tril(knnmat,-1), directed=None)
print ('Nodes=' + str(G.vcount()) + 'Edges=' + str(G.ecount()))
if not (len(G.components().sizes())==1):
	print ('Multiple components in the knn' + str(len(G.components().sizes())))

print ('louvain clnum=' + str(cl))
clst = []
res = 0
(clst,res) = getNClusters(G,cl)
if len(clst)>0:
	cname = 'knn' + str(k) + '_cl' + str(cl)
	rmat.insert(len(rmat.columns),cname,clst)
	resT.loc[cname] = res
else:
	print ('Failed to reach to the desired number of clusters')

rmat.to_csv(rmat_name,sep='\t')
resT.to_csv(resT_name,sep='\t')


