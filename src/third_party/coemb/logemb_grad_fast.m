function [f,g] = logemb_grad_fast(x,params)

dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;

[i,j] = max(pxy0);

[phi,psi] = read_legrad_x(x,params);
phi_tr = phi';
psi_tr = psi';
px0 = sum(pxy0,2);
py0 = sum(pxy0,1);

if 0
tic
% Generate model distribution
d = my_dist(phi,psi_tr);

pxy = exp(-d);
Z = sum(pxy(:));
pxy = pxy/Z;

px = sum(pxy,2);
py = sum(pxy,1);


gphi1 = diag(px)*phi-pxy*psi - (diag(px0)*phi-pxy0*psi);
%gpsi2 = diag(py)*psi-pxy'*phi - (diag(py0)*psi-pxy0'*phi);
gpsi1 = psi_tr*diag(py)-phi_tr*pxy - (psi_tr*diag(py0)-phi_tr*pxy0);
gpsi1 = gpsi1';
%[gphi,gpsi] = logemb_grad_c(phi,psi,phi_tr,psi_tr,px0,py0,phi_exp_d,psi_exp_d);
%gpsi = gpsi';
    
g = [-gphi1(:);-gpsi1(:)];

f = log(Z)+sum(pxy0(:).*d(:));
fprintf('Mlab took %g\n',toc);
else

%tic
% Test c code
psi_exp_d = pxy0*psi;
phi_exp_d = phi_tr*pxy0;
[gphi,gpsi,logpxy] = logemb_grad_c(phi,psi,phi_tr,psi_tr,px0,py0,phi_exp_d,psi_exp_d);
gpsi = gpsi';
f = -sum(pxy0(:).*logpxy(:));
g = [-gphi(:);-gpsi(:)];
%fprintf('C took %g\n',toc);
end
%fprintf('Grad diff=%g Func Diff=%g\n',max(abs(g-g1)),abs(f-f1));

%11;