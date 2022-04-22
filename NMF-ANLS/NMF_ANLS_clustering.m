function NMF_ANLS_clustering(alldata,rownames,colnames,k,outdirname)
% Written by Junha Shin (junha.shin@wisc.edu)
% Copyright (c) 2021 by Junha Shin and Sushmita Roy,
% Wisconsin Institute for Discovery, University of Wisconsin-Madison
% All rights reserved.
%
% INPUT alldata : [cells x genes] matrix, matrix value only without row and column headers
% INPUT rownames : [cells x 1] cell array, usually cell (row) names
% INPUT colnames : [1 x genes] cell array, usually gene (column) names
% INPUT k : integer, k number for NMF and clusterings
% INPUT outdirname : string, output directory naming
% OUTPUT : [outdirname]_K##/ directory, contains NMF result matrices (nmf_outs.mat) and clustering results as text

rng('shuffle');
addpath('nmf_bpas');	% locate the nmf_bpas directory
genenames=colnames;	% 1D-row array
cellnames=rownames;	% 1D-column array

% output naming
outsubdir=sprintf('%s/k%02d',outdirname,k)
mkdir(outsubdir);

% NMF/ANLS
disp (' - Calculating sparse NMF')
[u,v]=nmf(alldata,k,'type','sparse','nnls_solver','as');	% NMF-ANLS

%save results as nmf_outs.mat
outs{1}=u; outs{2}=v;
save(sprintf('%s/nmf_outs',outsubdir),'outs');

% u matrix maximal loadings
[ig,uid]=max(u');
fid=fopen(sprintf('%s/cell_cluster_maxdim_k%02d.txt',outsubdir,k),'w');
for i=1:length(uid)
	fprintf(fid,'%s\t%d\n',cellnames{i},uid(i));
end
fclose(fid);

% u matrix factor normalization
nmat = u;
for i=1:size(nmat,1)
	n=norm(nmat(i,:));
	if n==0
		n=1;
	end
	nmat(i,:) = nmat(i,:)/n;
end

% k-means (same k as nmf) of u matrix
outkmeans=sprintf('%s/cell_cluster_kmeans_k%02d.txt',outsubdir,k)
cidx=kmeans(nmat,k,'replicates',100);
result_kmeans=table(cellnames,cidx);
writetable(result_kmeans,outkmeans,'Delimiter','\t','WriteVariableNames',false);

return
