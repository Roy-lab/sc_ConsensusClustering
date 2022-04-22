function simmat_calculation_for_louvainPy(input_mat_file,using_colnum,dist,outname)
% this will generate sparse format simularity (affinity) matrix with index file

% loading matrix info
matf = importdata(input_mat_file);
matrix = matf.data(:,1:using_colnum);
rownames = matf.textdata(2:end,1);
colnames = matf.textdata(1,2:end);
size(matrix)

% calculate distance / similarity
disp(' - distance/similarity calculation')
if strcmp(dist,'euc')
	dmat = squareform(pdist(matrix));
elseif strcmp(dist,'pcc')
	dmat = squareform(pdist(matrix,'correlation'));
	dmat(isnan(dmat)) = 1;
elseif strcmp(dist,'cos')
	dmat = squareform(pdist(matrix,'cosine'));
	dmat(isnan(dmat)) = 1;
else
	error('ERROR: use proper code for distance metric');
end
smat = exp(-(dmat).^2 ./ (2*(std2(dmat))^2));	% affinity matrix
size(smat)

% write out
%[row col v] = find(sparse(smat));	% formatting to sparse matrix
[row col v] = find(sparse(tril(smat)));	% take only tril of smat, formatting to sparse matrix
simtrilmat = [row col v];
%outfile = sprintf('simmat_%s_%s_sparsetrilmat.txt',outname,dist)
%dlmwrite(outfile,[row col v],'delimiter','\t');
outfile = sprintf('simmat_%s_%s_sparsetrilmat.mat',outname,dist)
save(outfile,'simtrilmat','-v7.3')
outfile = sprintf('simmat_%s_%s_sparsetrilmat_index.txt',outname,dist)
idx = 0:(length(rownames)-1);
rtab = table(rownames,idx');
writetable(rtab,outfile,'delimiter','\t','writerownames',false,'writevariablenames',false)

quit

