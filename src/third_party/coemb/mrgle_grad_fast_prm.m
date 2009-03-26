function [f,gphi,gpsi,ga,gb,logpxy] = mrgle_grad_fast_prm(phi,psi,a,b,pxy0,px0,py0,b_cond)

phi_tr = phi';
psi_tr = psi';

[i,j,v] = find(pxy0);
psi_exp_d = pxy0*psi;
phi_exp_d = phi_tr*pxy0;

%[gphi,gpsi,ga,gb,f] = ...
%    mway_code_lowspace(full(phi),full(psi),full(phi_tr),full(psi_tr),full(px0),full(py0),full(phi_exp_d),full(psi_exp_d'),full(a),full(b),b_cond,i,j,v);
tic
[gphi,gpsi,ga,gb,f] = code_grad(full(phi),full(psi),full(phi_tr),full(psi_tr),full(px0),full(py0),full(phi_exp_d),full(psi_exp_d'),full(a),full(b),b_cond,i,j,v);

f = -f;


gphi = -gphi';
gpsi = -gpsi';
ga = ga;
gb = gb;

