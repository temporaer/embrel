function [test_perp,unseen_perp,lam] = code_delinterp(phi,psi,a,b,train,heldout,test,unseen_test)

[iho,jho,vho] = find(heldout);

NX = size(phi,1);
NY = size(psi,1);
dim = size(phi,2);
phi_exp_d = zeros(NY,dim);
psi_exp_d = zeros(NX,dim);
px0 = zeros(NX,1);
py0 = zeros(NY,1);
b_cond = 0;

%[gphi,gpsi,logpxy,ga,gb] =  mrgle_grad_wcond_c(phi,psi,phi',psi',px0,py0,phi_exp_d,psi_exp_d,a,b,b_cond);
px0 = sum(heldout,2);
py0 = sum(heldout,1);
[i,j,v] = find(heldout);
[gphi,gpsi,ga,gb,f,ho_lik] = condcode_lowspace_getlik(phi,psi,phi',psi',px0,py0,phi_exp_d,psi_exp_d',a,b,1,i,j,v);


ho_lik = exp(ho_lik);

[lam,hist,perp] = delinterp(train,heldout,ho_lik',0);

px0 = sum(test,2);
py0 = sum(test,1);
[i,j,v] = find(test);
[gphi,gpsi,ga,gb,f,tst_lik] = ...
    condcode_lowspace_getlik(phi,psi,phi',psi',px0,py0,phi_exp_d,psi_exp_d',a,b,1,i,j,v);

tst_lik = exp(tst_lik);

[aa,bb,test_perp] = delinterp(train,test,tst_lik',lam);
test_perp = full(test_perp);

px0 = sum(unseen_test,2);
py0 = sum(unseen_test,1);
[i,j,v] = find(unseen_test);
[gphi,gpsi,ga,gb,f,unseen_lik] = ...
    condcode_lowspace_getlik(phi,psi,phi',psi',px0,py0,phi_exp_d,psi_exp_d',a,b,1,i,j,v);

unseen_lik = exp(unseen_lik);

[a,b,unseen_perp] = delinterp(train,unseen_test,unseen_lik',lam);

unseen_perp = full(unseen_perp);

%load train_mc5000                                   
%train = mc_dat;
%load validation_mc5000
%valid = mc_dat;
%[lam,hist] = code_delinterp(phi_w1,psi_w2,a,b,train,valid);