% loading data matrix
matf=importdata('example_input/example_matrix.txt');
mat=matf.data;
cellnames=matf.textdata(2:end,1);
genenames=matf.textdata(1,2:end);


% Step1: Generating consensus source clustering results
% As an example, here we generated k={3,4,5} with 1 random initialization
disp(' - Generating consensus source clustering results')
addpath('NMF-ANLS')
k_start=3; k_increase=1; k_end=5; randinit=1;
outdir='example_output';
outprefix=sprintf('%s/example',outdir);
source_NMF_generation(mat,cellnames,genenames,k_start,k_increase,k_end,randinit,outprefix)


% Step2: Generating consensus graph matrix
disp(' - Generating consensus graph matrix')
inkmeanslist=sprintf('%s/example_NMF_sources/list_kmeans_results.txt',outdir)
outgraphmat=sprintf('%s/example_consensus_graph_matrix.txt',outdir)
command=sprintf('./consensus_code/consensus_clusters_graph %s %s',inkmeanslist,outgraphmat);
system(command);


% Step3: Generating consensus NMF clusters
% As an example, here we generated k={3,4,5} consensus clusters
disp(' - Generating consensus NMF clusters')
matf=importdata(outgraphmat);
mat=matf.data;
cellnames=matf.textdata(2:end,1);
outresultdir=sprintf('%s/consensus_results',outdir);
for k=3:5
	NMF_ANLS_clustering(mat,cellnames,cellnames,k,outresultdir)
end


