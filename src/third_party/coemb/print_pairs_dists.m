function print_pairs_dists(filename,phi,genes,tab_mode)
%
%
%

if(nargin<4)
  tab_mode=0;
end


nx = length(phi);
d = pdist(phi,'euclidean');
fid=fopen(filename,'w');

if(tab_mode)
  % Print in matrix format 
  sqd = squareform(d);  
  % Print header row 
  fprintf(fid,'.');
  for(i_gene=1:nx)
    fprintf(fid,'\t%s',genes{i_gene});
  end  
  fprintf(fid,'\n');

  % Print content 
  for(i_gene=1:nx)
    fprintf(fid,'%s',genes{i_gene});    
    for(j_gene=1:nx)  
      if(j_gene<=i_gene)
	fprintf(fid,'\t0');	
      else
	fprintf(fid,'\t%7.4f', sqd(i_gene,j_gene));
      end
    end
    fprintf(fid,'\n');
  end
  fclose(fid)
else  
  % Print a list 
  i=0;
  for(i_gene=1:nx)
    for(j_gene=i_gene+1:nx)  
      i=i+1;
      fprintf(fid,'%s\t%s\t%7.4f\n',...
	      genes{i_gene},...
	      genes{j_gene},...    
	      d(i));
    end
  end
  fclose(fid)
end

return