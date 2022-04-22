# graph_louvain_clustering
- **required** : matlab, python
- **required python packages** : numpy, scipy, pandas, igraph, louvain, h5py

See the argument usages in **run_louvain_resolution.sh** for a single run. <br>
See the parameters and arguments in **iter_run_louvain_resolution.sh** for a iterative runs.

This scripts will calculate Louvain clustering by taking input data matrix. Any type of data matrix could be used and **rows** will be clustered. <br>
Briefly, **simmat_calculation_for_louvainPy.m** will generate temporary similarty matrix to be used in the following **louvain_from_sparse_simmat_by_resolution.py**, which will calculate Louvain clustering based on given k of knn and resolution. <br>
Note that the used louvain algorithm is based on python package named vtraag (https://github.com/vtraag/louvain-igraph) which is known to be performed best. <br>

