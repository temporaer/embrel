function [f,g] = logemb_grad(x,params)

dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;

[phi,psi] = read_legrad_x(x,params);

% Generate model distribution
d = dist(phi,psi');

pxy = exp(-d.^2);
pxy = pxy/sum(pxy(:));

px = sum(pxy,2);
py = sum(pxy,1);
px0 = sum(pxy0,2);
py0 = sum(pxy0,1);

gphi = diag(px)*phi-pxy*psi - (diag(px0)*phi-pxy0*psi);
gpsi = diag(py)*psi-pxy'*phi - (diag(py0)*psi-pxy0'*phi);

g = [-gphi(:);-gpsi(:)];

f = -sum(sum(pxy0.*log(pxy)));