function [phi,psi,lik] = mrgle_grad_dim(pxy0,dim_list,max_repeats,init)
%
%
%

%pxy0 = pxy0';

if(exist('max_repeats')==0)
  max_repeats=1;
end
t0=cputime;
[NX,NY] = size(pxy0);
fprintf('Run gradient ascent (marg logemb_grad_dim)\n');

for curr_dim = dim_list
  for(i_repeat=1:max_repeats)
    fprintf('i_repeat=%d   cputime=%4.2fmin\n',i_repeat,(cputime-t0)/60);    
    if nargin<4
      phi0 = rand(NX,curr_dim);
      psi0 = rand(NY,curr_dim);    
      %    psi0 = rand(NY,curr_dim)*10-10;
      a0 = rand(NX,1);
      b0 = rand(NY,1);
    else
        phi0 = init.phis;
	psi0 = init.psis;
	a0 = init.a;
	b0 = init.b;
    end	
      
    x0 = [phi0(:);psi0(:);a0(:);b0(:)];
    cg_params.NX = NX;
    cg_params.NY = NY;
    cg_params.dim = curr_dim;
    cg_params.pxy0 = pxy0;
    
    opts = -1*zeros(1,9);
    x = conj_grad('mrgle_grad_fast',cg_params,x0,opts);
%  x = bfgswopt(x0,'logemb_grad',cg_params,1e-6,1e5);
%    options = optimset('GradObj','on','Display','iter');
%    x = fminunc('logemb_grad',x0,options,cg_params);
    
    [phis{i_repeat},psis{i_repeat}] = read_mrgle_x(x,cg_params);
    liks(i_repeat) = -full(mrgle_grad_fast(x,cg_params));
    fprintf('Lik=%g',liks(i_repeat));
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
