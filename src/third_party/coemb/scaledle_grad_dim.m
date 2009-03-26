function [phi,psi,lik,liks] = scaledle_grad_dim(pxy0,dim_list,beta,init)
%
%
%

%pxy0 = pxy0';

max_repeats=20;
t0=cputime;
[NX,NY] = size(pxy0);
fprintf('Run gradient ascent (logemb_grad_dim)\n');

cg_params.beta = beta;

for curr_dim = dim_list
  for(i_repeat=1:max_repeats)
    fprintf('i_repeat=%d   cputime=%4.2fmin\n',i_repeat,(cputime-t0)/60);    
    
    if nargin>=4
      phi0 = init.phis;
      psi0 = init.psis;
    else
      phi0 = 1e-4*rand(NX,curr_dim);
      psi0 = 1e-4*rand(NY,curr_dim);
    end
    a0 = log(sum(pxy0,2));rand(NX,1);
    b0 = log(sum(pxy0,1)); rand(NY,1);
    %        x0 = [xv0(:);yv0(:)];

    cg_params.NX = NX;
    cg_params.NY = NY;
    cg_params.dim = curr_dim;
    cg_params.px0 = sum(pxy0,2);
    cg_params.py0 = sum(pxy0,1);    
    cg_params.oprod = cg_params.px0*cg_params.py0;
    cg_params.lik0 = sum(pxy0(:).*log(cg_params.oprod(:)));
    cg_params.pxy0 = pxy0;
    cg_params.beta = beta;
%    cg_params.nns = nns';

    opts = -1*zeros(1,9);

%    x0 = [phi0(:);psi0(:);beta];
    x0 = [phi0(:);psi0(:)];
     
     x = conj_grad('scaledle_grad',cg_params,x0,opts);
%     x = proj_grad('scaledle_grad','',cg_params,x0,5000);
%    [phis{i_repeat},psis{i_repeat}] =
%    read_scaledle_x(x,cg_params);
    [phis{i_repeat},psis{i_repeat}] = read_legrad_x(x,cg_params);
%    [phis{i_repeat},psis{i_repeat}] = read_scaledle_x(x,cg_params);
    liks(i_repeat) = -scaledle_grad(x,cg_params);    
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
