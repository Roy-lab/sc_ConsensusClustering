function source_NMF_generation(matrix,rownames,colnames,k_start,k_inc,k_end,I_max,outprefix)
% matrix : [cells x genes] matrix, matrix value only without row and column headers
% rownames : [cells x 1] cell array, usually cell (row) names
% colnames : [1 x genes] cell array, usually gene (column) names
% k_start : integer, starting k number for the iteration
% k_inc : integer, incremental number for k for the iteration 
% k_end : integer, ending k number for the iteration
% I_max : integer, number of random initialization
% outprefix : string, outname prefix, result directory name will be "[XX]_NMF_sources/"

addpath('NMF-ANLS/nmf_bpas/');	% NMF-ANLS code

list_maxdim={};
list_kmeans={};
for I=1:I_max
for k=k_start:k_inc:k_end       % e.g. for k=20:10:30
        outsubdir=sprintf('%s_NMF_sources/k%02d-I%02d',outprefix,k,I)
        mkdir(outsubdir);

        % NMF-ANLS
	rng('shuffle');
        [u,v]=nmf(matrix,k,'type','sparse','nnls_solver','as');        % row=cells, col=genes

        %save allmatrix as .mat
        outs{1}=u; outs{2}=v;
        save(sprintf('%s/nmf_outs',outsubdir),'outs');

        % maxdim
        [ig,uid]=max(u');
        fid=fopen(sprintf('%s/cellclust_maxdim.txt',outsubdir),'w');
        for i=1:length(uid)
                fprintf(fid,'%s\t%d\n',rownames{i},uid(i));
        end
        fclose(fid);
	list_maxdim{end+1}=sprintf('%s/cellclust_maxdim.txt',outsubdir);

        % k-means with factor normalization
        umat = u;
        for i=1:size(umat,1)
                n=norm(umat(i,:));
                if n==0, n=1; end
                umat(i,:) = umat(i,:)/n;
        end
        kk=k;
        outkmeans=sprintf('%s/cellclust_kmeans.txt',outsubdir);
        cidx=kmeans(umat,kk,'replicates',100);
        result_kmeans=table(rownames,cidx);
        writetable(result_kmeans,outkmeans,'Delimiter','\t','WriteVariableNames',false);
	list_kmeans{end+1}=sprintf('%s/cellclust_kmeans.txt',outsubdir);
end
end

outfilename=sprintf('%s_NMF_sources/list_maxdim_results.txt',outprefix)
rtab=table(list_maxdim');
writetable(rtab,outfilename,'WriteRowNames',false,'WriteVariableNames',false);

outfilename=sprintf('%s_NMF_sources/list_kmeans_results.txt',outprefix)
rtab=table(list_kmeans');
writetable(rtab,outfilename,'WriteRowNames',false,'WriteVariableNames',false);

return

