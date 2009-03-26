function [f,g] = copula_logemb_grad_fast(x,params)

dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;
px0 = params.px0;
py0 = params.py0;

[phi,psi] = read_legrad_x(x,params);
phi_tr = phi';
psi_tr = psi';


%tic
% Test c code
psi_exp_d = pxy0*psi;
phi_exp_d = phi_tr*pxy0;

[gphi,gpsi,logpxy] = ...
    copula_logemb_grad_c(phi,psi,phi_tr,psi_tr,px0,py0,phi_exp_d,psi_exp_d);
gpsi = gpsi';
f = -sum(pxy0(:).*logpxy(:))+params.lik0;
g = [-gphi(:);-gpsi(:)];
