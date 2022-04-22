set -u
# required : matlab, python
# required python packages : numpy, scipy, pandas, h5py, igraph, louvain

matfile=$1	# e.g. vv_submat.txt
setname=$2	# e.g. vv
knnk=$3		# e.g. 20
RES=$4		# e.g. 0.1
USING_COLNUM=$5	# Number of data dimensions (i.e. columns) for calculating Euc dist. of
                # similarity matrix. This distance will be the edge weights for taking
                # knn graph. Use full count of columns for the case full-matrix
                # (rather than dim reduced matrix) for the graph generation. See
                # simmat_calculation_for_louvainPy.m in details.

matlab -r simmat_calculation_for_louvainPy\(\'$matfile\',$USING_COLNUM,\'euc\',\'$setname\'\)
	
outspsimmat=simmat_${setname}_euc_sparsetrilmat.mat
outspsimidx=simmat_${setname}_euc_sparsetrilmat_index.txt

python louvain_from_sparse_simmat_by_resolution.py $outspsimmat $outspsimidx $knnk $RES $setname

rm $outspsimmat $outspsimidx

