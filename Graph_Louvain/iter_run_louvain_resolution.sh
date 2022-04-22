set -u
# required : matlab, python
# required python packages : numpy, scipy, pandas, h5py, igraph, louvain

inputmatrix=$1  # e.g. NMF u matrix or full data matrix
outprefix=$2    # output name prefix, e.g. "vv_nmfk50"

knnks=(10 20 30 40 50)			# list of k numbers for knn graph
resolutions=(0.1 0.3 0.5 0.7 1.0)	# list of Louvain resolutions
using_colnum=50	# Number of data dimensions (i.e. columns) for calculating Euc dist. of 
		# similarity matrix. This distance will be the edge weights for taking 
		# knn graph. Use full count of columns for the case full-matrix 
		# (rather than dim reduced matrix) for the graph generation. See 
		# simmat_calculation_for_louvainPy.m in details.


for knnk in "${knnks[@]}"
do
resultdir=${outprefix}_louvain_withRes
mkdir -p ${resultdir}/knn$knnk/

matlab -r simmat_calculation_for_louvainPy\(\'$inputmatrix\',$using_colnum,\'euc\',\'$outprefix\'\)
outspsimmat=simmat_${outprefix}_euc_sparsetrilmat.mat
outspsimidx=simmat_${outprefix}_euc_sparsetrilmat_index.txt

	for res in "${resolutions[@]}"
	do	
	mkdir ${resultdir}/knn$knnk/res$res/

	python louvain_from_sparse_simmat_by_resolution.py $outspsimmat $outspsimidx $knnk $res $outprefix
	outfilename=${outprefix}_louvain_knn${knnk}_res${res}.txt
	# v2a_nmfk50_louvain_knn10_res0.1.txt

	mv $outfilename ${outprefix}_louvain.txt
	# v2a_nmfk50_louvain.txt
	mv ${outprefix}_louvain.txt ${resultdir}/knn$knnk/res$res/
	done

rm $outspsimmat $outspsimidx
done

