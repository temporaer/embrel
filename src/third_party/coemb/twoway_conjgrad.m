function [phi,psis,as,bs,lik]= twoway_conjgrad(pxy0s,dim_list,options)

[NX,NY] = size(pxy0s{1});
% Initailize phi0 to something
options.phi_init = rand(NX,dim_list);
options.psi_init{1} = rand(NY,dim_list);
options.a_init = rand(NX,1);
options.b_init{1} = rand(NY,1);

options.n_restarts = 1;
options.b_keep_phi = 1;
options.b_keep_psi = 0;

for i=1:100
  [phi,psis,as,bs,lik] = code(pxy0s,dim_list,options);
  options.psi_init{1} = psis{dim_list};
  options.phi_init=phi{dim_list};
  options.a_init = as;
  options.b_init = bs;
  options.b_keep_phi = 0;
  options.b_keep_psi = 1;
  
  [phi,psis,as,bs,lik] = code(pxy0s,dim_list,options);  
  options.phi_init = phi{dim_list};
  options.psi_init{1} = psis{dim_list};
  options.b_keep_phi = 1;
  options.b_keep_psi = 0;
  options.a_init = as;
  options.b_init = bs;
  
  fprintf('Iter=%d Lik=%g\n',i,lik);
end



