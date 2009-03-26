function m = vec_to_smat_fast(x,N)

m = vec_to_smat_c(full(x),N);
m = m + tril(m,-1)';
if 0
m = zeros(N,N);
xind = 1;

for i=N:-1:1
  rowind = N-i+1;
  m(rowind,rowind:end) = x(xind:xind+i-1)';
  xind = xind+i;
end

m = m + triu(m,1)';
end