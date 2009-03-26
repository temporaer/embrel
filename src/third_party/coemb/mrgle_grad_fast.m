function [f,g] = mrgle_grad_fast(x,params)

dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;

[phi,psi,a,b] = read_mrgle_x(x,params);
phi_tr = phi';
psi_tr = psi';
px0 = sum(pxy0,2);
py0 = sum(pxy0,1);


%tic
% Test c code
psi_exp_d = pxy0*psi;
phi_exp_d = phi_tr*pxy0;
[gphi,gpsi,logpxy,ga,gb] = mrgle_grad_c(phi,psi,phi_tr,psi_tr,px0,py0,phi_exp_d,psi_exp_d,a,b);

gpsi = gpsi';
f = -sum(pxy0(:).*logpxy(:));
%gb = gb*0;
%ga = ga*0;
g = [-gphi(:);-gpsi(:);ga(:);gb(:)];
%fprintf('C took %g\n',toc);

%fprintf('Grad diff=%g Func Diff=%g\n',max(abs(g-g1)),abs(f-f1));

%11;
