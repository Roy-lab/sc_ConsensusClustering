
**NMF_ANLS_clustering.m**
----------
**NOTE: this script works doing clustering with 1 specific k number.**<br>
**usage: NMF_ANLS_clustering(alldata,rownames,colnames,k,outdirname)**
* alldata : [cells x genes] matrix, matrix value only without row and column headers
* rownames : [cells x 1] cell array, usually cell (row) names
* colnames : [1 x genes] cell array, usually gene (column) names
* k : integer, k number for NMF and clusterings
* outdirname : string, output directory naming, result will be [outdirname]/k## directory, contains NMF result matrices (nmf_outs.mat) and clustering results as text

**source_NMF_generation.m**
----------
**NOTE: this script basically works doing clustering with multiple times.**<br>
**usage: source_NMF_generation(matrix,rownames,colnames,k_start,k_inc,k_end,I_max,outprefix)**
* matrix : [cells x genes] matrix, matrix value only without row and column headers
* rownames : [cells x 1] cell array, usually cell (row) names
* colnames : [1 x genes] cell array, usually gene (column) names
* k_start : integer, starting k number for the iteration
* k_inc : integer, incremental number for k for the iteration 
* k_end : integer, ending k number for the iteration
* I_max : integer, number of random initialization
* outprefix : string, outname prefix, result directory name will be "[XX]_NMF_sources/"
