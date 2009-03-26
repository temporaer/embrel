function loc = find_incell(str,cellarr,multi)

loc=[];
if nargin<3
  multi = 0;
end

for i=1:length(cellarr)
    if strcmp(cellarr{i},str)
      if ~multi
	loc = i;
        return;
      end
      loc(end+1) = i;
    end
end
if ~multi
  loc=-1;
end
