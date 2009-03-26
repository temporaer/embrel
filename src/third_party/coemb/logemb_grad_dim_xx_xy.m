function [phi,psi,lik] = logemb_grad_dim_xx_xy(pxy0,pxx0,dim_list,max_repeats)
%
%
%

%pxy0 = pxy0';

if(exist('max_repeats')==0)
  max_repeats=10;
end
t0=cputime;
[NX,NY] = size(pxy0);
fprintf('Run gradient ascent (logemb_grad_dim)\n');

for curr_dim = dim_list
  for(i_repeat=1:max_repeats)
    fprintf('i_repeat=%d   cputime=%4.2fmin\n',i_repeat,(cputime-t0)/60);    
    xv0 = 1e-3*rand(NX,curr_dim);
    yv0 = 1e-3*rand(NY,curr_dim);
    
    x0 = [xv0(:);yv0(:)];
    cg_params.NX = NX;
    cg_params.NY = NY;
    cg_params.dim = curr_dim;
    cg_params.pxy0 = pxy0;
    cg_params.pxx0 = pxx0;
    cg_params.lambda = 1;
    
    opts = -1*zeros(1,9);
    x = conj_grad('logemb_grad_xx_xy',cg_params,x0,opts);
    
    [phis{i_repeat},psis{i_repeat}] = read_legrad_x(x,cg_params);
    liks(i_repeat) = -logemb_grad_xx_xy(x,cg_params);
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

%tmp_phi = phi;
%tmp_psi = psi;
%phi = tmp_psi;
%psi = tmp_phi;

return
