function [f,g] = scaledle_grad(x,params)

dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;
beta = params.beta;
oprod = params.oprod;

%[v,j] = max(pxy0(:));
%[i,j] = ind2sub(size(pxy0),j);

%[phi,psi,beta] = read_scaledle_x(x,params);
%beta = 1e4;
[phi,psi] = read_legrad_x(x,params);
phi_tr = phi';
psi_tr = psi';

% Generate model distribution
d = my_dist(phi,psi_tr);
dsum = sum(d(:));

phi = phi/dsum;
psi = psi/dsum;
phi_tr = phi_tr/dsum;
psi_tr = psi_tr/dsum;

%pxy = oprod.*exp(-beta*d/dsum);
pxy = exp(-beta*d/dsum);
%mx = max(d(:)/dsum)*beta;
%logpxy = mx-beta*d/dsum-log(sum(sum(exp(-d/dsum*beta+mx))));%logpxy = mx-beta*d/dsum-log(sum(sum(oprod.*exp(-d/dsum*beta+mx))))+log(oprod);

Z = sum(pxy(:));
pxy = pxy/Z;
%pxy = exp(logpxy);

px = sum(pxy,2);
py = sum(pxy,1);
px0 = sum(pxy0,2);
py0 = sum(pxy0,1);

sumpsi = sum(psi,1);
sumphi = sum(phi,1);

gphi = diag(px)*phi-pxy*psi - (diag(px0)*phi-pxy0*psi);
gphi = gphi +(NY*phi-repmat(sumpsi,NX,1))*sum(pxy0(:).*d(:)/dsum);
gphi = gphi -(NY*phi-repmat(sumpsi,NX,1))*sum(pxy(:).*d(:)/dsum);
%gpsi2 = diag(py)*psi-pxy'*phi - (diag(py0)*psi-pxy0'*phi);
gpsi = psi_tr*diag(py)-phi_tr*pxy - (psi_tr*diag(py0)-phi_tr*pxy0);
gpsi = gpsi' +(NX*psi-repmat(sumphi,NY,1))*sum(pxy0(:).*d(:)/dsum);
gpsi = gpsi  -(NX*psi-repmat(sumphi,NY,1))*sum(pxy(:).*d(:)/dsum);

%gbeta = 1*sum((pxy(:)-pxy0(:)).*d(:)/dsum);
%gbeta = 0;
%g = [-gphi(:);-gpsi(:);-gbeta];
g = [-gphi(:);-gpsi(:)];

%f = -sum(pxy0(:).*logpxy(:))+params.lik0;
f = -sum(pxy0(:).*log(pxy(:)+eps));

