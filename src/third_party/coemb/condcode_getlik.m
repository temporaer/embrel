function p = condcode_getlik(phi,psi,a,b,dat)

NX = size(phi,1);
NY = size(psi,1);
dim = size(phi,2);
phi_exp_d = zeros(NY,dim);
psi_exp_d = zeros(NX,dim);
px0 = zeros(NX,1);
py0 = zeros(NY,1);
b_cond = 0;

%[gphi,gpsi,logpxy,ga,gb] =  mrgle_grad_wcond_c(phi,psi,phi',psi',px0,py0,phi_exp_d,psi_exp_d,a,b,b_cond);
px0 = sum(dat,2);
py0 = sum(dat,1);
[i,j,v] = find(dat);

[gphi,gpsi,ga,gb,f,lik] = condcode_lowspace_getlik(phi,psi,phi',psi',px0,py0,phi_exp_d,psi_exp_d',a,b,1,i,j,v);

lik = exp(lik);

p = full(sparse(i,j,lik));