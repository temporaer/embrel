function [phi,psi,lik,liks] = lenn_grad_dim(pxy0,dim_list,init)
%
%
%

pxy0 = pxy0';

if(exist('max_repeats')==0)
  max_repeats=1;
end
t0=cputime;
[NX,NY] = size(pxy0);
fprintf('Run gradient ascent (logemb_grad_dim)\n');

if exist('beta','var')
    cg_params.beta = beta;
else
    cg_params.beta = 1;
end

for curr_dim = dim_list
  for(i_repeat=1:max_repeats)
    fprintf('i_repeat=%d   cputime=%4.2fmin\n',i_repeat,(cputime-t0)/60);    
    
    if nargin<3
      phi0 = 1e-3*rand(NX,curr_dim);
      psi0 = 1e-3*rand(NY,curr_dim);
      a0 = log(sum(pxy0,2));rand(NX,1);
      b0 = log(sum(pxy0,1)); rand(NY,1);
      %        x0 = [xv0(:);yv0(:)];
    else
      phi0 = init.phis;
      psi0 = init.psis;
      a0 = init.a0;
      b0 = init.b0;
    end	
    cg_params.NX = NX;
    cg_params.NY = NY;
    cg_params.dim = curr_dim;
    cg_params.x0 = sum(pxy0,2);
    cg_params.y0 = sum(pxy0,1);    
    cg_params.oprod = cg_params.x0*cg_params.y0;
    cg_params.lik0 = sum(pxy0(:).*log(cg_params.oprod(:)));
    cg_params.pxy0 = pxy0;
%    cg_params.nns = nns';
    
    opts = -1*zeros(1,9);
%    x = conj_grad('logemb_grad_fast',cg_params,x0,opts);
%    x = conj_grad('copula_logemb_grad',cg_params,x0,opts);
%     phi0 = 1e-5*rand(NX,curr_dim);
%     psi0 = 1e-5*rand(NY,curr_dim);
     
     
%     x0 = [phi0(:);psi0(:);a0(:);b0(:)];
     x0 = [phi0(:);psi0(:)];

%    x = conj_grad('mrgle_grad_fast',cg_params,x0,opts);
    x = conj_grad('logemb_grad_fast',cg_params,x0,opts);
%    cg_params.sig  = sig;
%    cg_params.bet  = bet;    
%    x = conj_grad('prmod_logemb_grad',cg_params,x0,opts);
%    x = bfgswopt(x0,'lenn_grad',cg_params,1e-6,100);
%    x = conj_grad('lenn_grad',cg_params,x0,opts);    
%   x = proj_grad('copula_logemb_grad','logemb_proj',cg_params,x0,10000);    
%    x = proj_conj_grad('copula_logemb_grad','logemb_proj',cg_params,x0,opts);    
%    x = proj_grad('copula_logemb_grad','',cg_params,x0,1e5);    
%  x = bfgswopt(x0,'lenn_grad',cg_params,1e-6,1e5);
%    options = optimset('GradObj','on','Display','iter');
%    x = fminunc('logemb_grad',x0,options,cg_params);
    
    [phis{i_repeat},psis{i_repeat}] = read_legrad_x(x,cg_params);
%     [phis{i_repeat},psis{i_repeat},a,b] = ...
%	 read_mrgle_x(x,cg_params);

    liks(i_repeat) = -logemb_grad_fast(x,cg_params);    
%    liks(i_repeat) = -copula_logemb_grad(x,cg_params);
%    liks(i_repeat) = -prmod_logemb_grad(x,cg_params);    
%    liks(i_repeat) = -mrgle_grad_fast(x,cg_params);    
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

tmp_phi = phi;
tmp_psi = psi;
phi = tmp_psi;
psi = tmp_phi;

return
