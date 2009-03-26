d = dir('*.c');
for di=1:length(d)
  eval(sprintf('mex %s',d(di).name));
end
