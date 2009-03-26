function d = kappav(p)

n = length(p)+1;

diagp = diag(p);
e = ones(n-1,1);

x = -1/(n+sqrt(n)); y = -1/sqrt(n); 
% Calculate distances between new point and all the rest
const = (x-y)^2*sum(p(:));

newd = diagp+2*(x-y)*sum(p,2)+const;

d = zeros(n,n);
% Creates distances for all n-1 points
d(2:end,2:end) = diagp*e'+e*diagp'-2*p;
d(1,2:end) = newd';
d(2:end,1) = newd;

