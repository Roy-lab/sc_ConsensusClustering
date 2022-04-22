import sys, re, os
import pandas as pd
import scipy as sp
import scipy.sparse
import numpy as np
import h5py
import igraph as ig
import louvain

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


## load parameters
sparse_simmat_file = sys.argv[1]
sparse_simmat_idx = sys.argv[2]
knn_k = int(sys.argv[3])
louvain_resolution = float(sys.argv[4])
outname = sys.argv[5]

## load data
#spsimmat = np.genfromtxt(sparse_simmat_file,delimiter='\t')
arrays = {}
f = h5py.File(sparse_simmat_file,'r');
for k, v in f.items():
	arrays[k] = np.array(v)

spsimmat = arrays['simtrilmat'].T
spsimmatidx = np.genfromtxt(sparse_simmat_idx,dtype=str,delimiter='\t')
row = spsimmat[:,0] -1
row = row.astype(int)
col = spsimmat[:,1] -1
col = col.astype(int)
val = spsimmat[:,2]
simmat = sp.sparse.coo_matrix((val, (row, col)), shape=(len(spsimmatidx),len(spsimmatidx))).toarray()
simmat = simmat + simmat.T - np.diag(simmat.diagonal())	# symmetricize/restore the trilmat

## generate knn graph
k = knn_k
print ('kNN= ' + str(k))
knn_binmat = knn_bin_adjmat(simmat,k)

## assign edge weights
simdistmat = 1-simmat
knnmat=np.multiply(simdistmat,knn_binmat)# similarity=1 is assigned as "0" in this matrix
simdistmin=np.min(simdistmat[np.nonzero(simdistmat)])   # non-zero min distance
simdist_smallval = simdistmin * 0.1     # multiplying 0.1 to assign most small distance
print ('smallest distance=' + str(simdist_smallval))
smallvalidx=np.multiply((simdistmat==0),(knn_binmat==1))# index of similarity=1 edges
smallvalmat=np.multiply(simdist_smallval,smallvalidx)   # matrix of similarity=1 edges
knnmat = knnmat + smallvalmat

## generate graph
G = get_igraph_from_adjacency(np.tril(knnmat,-1), directed=None)
print ('Nodes= ' + str(G.vcount()) + ', Edges= ' + str(G.ecount()))
if not (len(G.components().sizes())==1):
	print ('Multiple components in the knn: ' + str(len(G.components().sizes())))
	#print ('Skipping knn= ' + str(k))
	#sys.exit()

## louvain clustering
clst = []
partition_type = louvain.RBConfigurationVertexPartition
partition_kwargs = dict()
partition_kwargs["resolution_parameter"] = louvain_resolution
weights = np.array(G.es["weight"]).astype(np.float64)
partition_kwargs["weights"] = weights
part = louvain.find_partition(G,partition_type,**partition_kwargs)
clusters = np.array(part.membership)
this_clusters = len(np.unique(clusters))
print('got ' + str(this_clusters) + ' cluster(s) at resolution ' + str(louvain_resolution))
clst = clusters+1

## write out
result = pd.DataFrame(data=clst,index=spsimmatidx[:,0],columns=[('louvain_' + str(louvain_resolution))])
outfile = outname + '_louvain_knn' + str(knn_k) + '_res' + str(louvain_resolution) + '.txt'
result.to_csv(outfile,sep='\t',header=0)


