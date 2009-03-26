function op = get_kappa_op_inds(n,xi,yi,ws)

% Find how Kappa operates linearily to transform a PSD (n-1) matrix
% into a EDM
% Do this by running Kappa on the standard basis of p
A=[];


lininds = sub2ind([n,n],xi,yi);

ind = 0;

x = -1/(n+sqrt(n));
y = -1/sqrt(n);
a_diff = 2*(x-y)^2;
b_diff = 2*(x+1-y)*(x-y);
a_eq = (x-y)^2;
b_eq = (x+1-y)^2;

op = sparse(length(lininds),sum(1:n-1));

% Only go over half the matrix
for i=1:n-1
  i  
  for j=i:n-1
    ind = ind+1;
    d = sparse(n,n);
    if i~=j
      d(i+1,j+1)=-2;
      d(1,2:end)=a_diff;
      d(1,i+1) = b_diff;
      d(1,j+1) = b_diff;
    else
     d(1,2:end) = a_eq;
      d(1,i+1) = b_eq;
      d(2:i,j+1) = 1;
      d(j+1,j+2:end) = 1;
    end      
    op(:,ind) = d(lininds);
   end
end

NX = max(xi);
NY = max(yi)-NX;
eval(sprintf('save A_op%d_NX%d_NY%d A',n,NX,NY));
