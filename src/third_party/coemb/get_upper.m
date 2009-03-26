function v =get_upper(m)

% Return upper triangular matrix, spread out as a vector

N = size(m,1);
xind = 1;
v = sparse(1,(N-1)*N/2);

for i=N-1:-1:1
  rowind = N-i;
  v(xind:xind+i-1) = m(rowind,rowind+1:end);
  xind = xind+i;
end
