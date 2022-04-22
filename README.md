# sc_ConsensusClustering
Scripts for unbiased clustering and consensus clustering

**NMF-ANLS/**
----------
matlab scripts for doing NMF-ANLS clustering.
Original nmf_bpas code was downloaded from: http://www.cc.gatech.edu/~hpark/software/nmf_bpas.zip

**Graph_Louvain/**
----------
python scripts for doing graph Louvain clustering.

**consensus_code/**
----------
A C++ code for generating consensus graph matrix of multiple clustering solutions.

**example_consensus_clustering_wrapper.m**
----------
**usage: matlab -r example_consensus_clustering_wrapper**<br>
An example wrapper script (matlab) for running through the procedures of the consensus clustering.
Here we assume that (1) generating sources based on NMF-ANLS kmeans and (2) generating results as NMF-ANLS.
