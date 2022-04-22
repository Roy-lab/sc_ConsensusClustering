# graph_louvain_clustering
* **required** : python
* **required python packages** : numpy, scipy, pandas, igraph, louvain
<br>

**data_matrix_louvain_clustering_by_k.py**
----------
* usage: python data_matrix_louvain_clustering_by_k.py [data_matrix] [k_of_knn] [target_number_of_clusters] [output_name]

**data_matrix_louvain_clustering_by_resolution.py**
----------
* usage: python data_matrix_louvain_clustering_by_k.py [data_matrix] [k_of_knn] [resolution_of_louvain] [output_name]
<br>

This scripts will calculate Louvain clustering by taking **input data matrix**.<br>
Any type of data matrix could be used. **Rows** will be clustered. <br>
Note that the used louvain algorithm is based on python package implemented by V. Traag (https://github.com/vtraag/louvain-igraph) which is known to be performed best. <br>
