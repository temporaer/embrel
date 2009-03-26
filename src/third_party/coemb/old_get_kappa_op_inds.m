function [grad,op] = get_kappa_op_inds(n,xi,yi,ws)

% Find how Kappa operates linearily to transform a PSD (n-1) matrix
% into a EDM
% Do this by running Kappa on the standard basis of p
A=[];
%grad = sparse(1,(n-1)^2);

lininds = sub2ind([n,n],xi,yi);

ind = 0;

% Only go over half the matrix
for i=1:n-1
  for j=i:n-1
    ind = ind+1;
    p=sparse(n-1,n-1);
    p(i,j) = 1;
    p(j,i) = 1;
    d = kappav(p);
    grad(ind) = dot(d(lininds),ws);
    op(:,ind) = d(lininds);
   end
end
