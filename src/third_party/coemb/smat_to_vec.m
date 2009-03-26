function v = smat_to_vec(m)

N = size(m,1);
xind = 1;
v = sparse(1,N*(N+1)/2);

for i=N:-1:1
  rowind = N-i+1;
  v(xind:xind+i-1) = m(rowind,rowind:end);
  xind = xind+i;
end
