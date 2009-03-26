function [phi,psi,lik,liks] = condcode_grad_dim(pxy0,dim_list,max_repeats,beta)
%
%
%

if nargin<3
  max_repeats=1;
end
t0=cputime;
[NX,NY] = size(pxy0);
fprintf('Run gradient ascent (logemb_grad_dim)\n');

if nargin>=4
    cg_params.beta = beta;
else
    cg_params.beta = 1;
end

for curr_dim = dim_list
  for(i_repeat=1:max_repeats)
    fprintf('i_repeat=%d   cputime=%4.2fmin\n',i_repeat,(cputime-t0)/60);    
    xv0 = rand(NX,curr_dim);
    yv0 = rand(NY,curr_dim);
    
    x0 = [xv0(:);yv0(:)];
    cg_params.NX = NX;
    cg_params.NY = NY;
    cg_params.dim = curr_dim;
    cg_params.pxy0 = pxy0;
    cg_params.px0 = sum(pxy0,2);
    cg_params.py0 = sum(pxy0,1);
    
    opts = -1*zeros(1,9);
%    x = conj_grad('logemb_grad_fast',cg_params,x0,opts);
    x = conj_grad('condcode_grad',cg_params,x0,opts);
%    x = conj_grad('copula_logemb_grad',cg_params,x0,opts);    
%  x = bfgswopt(x0,'logemb_grad',cg_params,1e-6,1e5);
%    options = optimset('GradObj','on','Display','iter');
%    x = fminunc('logemb_grad',x0,options,cg_params);
    
    [phis{i_repeat},psis{i_repeat}] = read_legrad_x(x,cg_params);
    liks(i_repeat) = -condcode_grad(x,cg_params);    
%    liks(i_repeat) = -copula_logemb_grad(x,cg_params);
    fprintf('Lik=%g\n',liks(i_repeat));
  end
  liks = full(liks)
  [max_lik,bst_repeat] = max(liks);
  bst_repeat = bst_repeat(1);
  fprintf('choose best score out of repeats\n');  
  fprintf('max_lik=%f bst_repeat=%d\n',max_lik,bst_repeat);    
  
  
  lik(curr_dim) = max_lik;
  phi = phis{bst_repeat};
  psi = psis{bst_repeat};  
end

return
