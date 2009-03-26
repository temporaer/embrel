function [phis,psi,as,bs] = read_mway_params_y(x,cg_params)

dim = cg_params.dim;
x_i = 1;

% Read the psi, which is common to all distributions
if ~cg_params.b_keep_psi
  sz = cg_params.NY*dim;
  psi = reshape(x(x_i:x_i+sz-1),cg_params.NY,dim);
  x_i = x_i+sz;
else
  psi = cg_params.psi0;
end

for i=1:length(cg_params.pxy0s)
  if ~cg_params.b_keep_phi
    sz = cg_params.NXs(i)*dim;  
    phis{i} = reshape(x(x_i:x_i+sz-1),cg_params.NXs(i),dim);  
    x_i = x_i + sz;
  else
    phis{i} = cg_params.psi0{i};
  end
  
  if ~cg_params.b_keep_a
      sz = cg_params.NXs(i);
      as{i} = x(x_i:x_i+sz-1);
      x_i = x_i+sz;
  else
      as{i} = cg_params.a0{i};
  end
  if ~cg_params.b_keep_b
      sz = cg_params.NY;  
      bs{i} = x(x_i:x_i+sz-1);
  else
      bs{i} = cg_params.b0{i};
  end
end
