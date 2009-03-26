function x = read_format_file(fn,pref)

fid = fopen(fn);
x=[];
while ~feof(fid)
  l = fgetl(fid);
  i = strfind(l,pref);
  if isempty(i)
    continue;
  end
  x(end+1) = str2num(l(i+length(pref):end));
end
fclose(fid);
