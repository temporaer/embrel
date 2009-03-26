function print_pairs_dists_to_tab(filename,phi,genes)
%
%
%

nx = length(phi);

d = pdist(phi,'euclidean');
sqd = squareform(d);
fid=fopen(filename,'w');
i=0;

fprintf(fid,'.');
for(i_gene=1:nx)
  fprintf(fid,'\t%s',genes{i_gene});
end

for(i_gene=1:nx)
  for(j_gene=i_gene+1:nx)  
    i=i+1;
    fprintf(fid,'%s\t%s\t%8.5f\n',...
	    genes{i_gene},...
	    genes{j_gene},...    
	    d(i));
  end
end
fclose(fid)

return