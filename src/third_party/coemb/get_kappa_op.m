function A = get_kappa_op(n)

% Find how Kappa operates linearily to transform a PSD (n-1) matrix
% into a EDM
% Do this by running Kappa on the standard basis of p
A=[];

for i=1:n-1
  for j=1:n-1
    p=sparse(n-1,n-1);
    p(i,j) = 1;
    d = kappav(p);
    A(:,end+1) = d(:);
    size(A,2)
  end
end
