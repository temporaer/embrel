function [phi,psis,as,bs] = read_mway_params(x,cg_params)

dim = cg_params.dim;
x_i = 1;

% Read the phi, which is common to all distributions
if ~cg_params.b_keep_phi
  sz = cg_params.NX*dim;
  phi = reshape(x(x_i:x_i+sz-1),cg_params.NX,dim);
  x_i = x_i+sz;
else
  phi = cg_params.phi0;
end

for i=1:length(cg_params.pxy0s)
  if ~cg_params.b_keep_psi
    sz = cg_params.NYs(i)*dim;  
    psis{i} = reshape(x(x_i:x_i+sz-1),cg_params.NYs(i),dim);  
    x_i = x_i + sz;
  else
    psis{i} = cg_params.psi0{i};
  end
  if ~cg_params.b_keep_a
      sz = cg_params.NX;
      as{i} = x(x_i:x_i+sz-1);
      x_i = x_i+sz;
  else
      as{i} = cg_params.a0{i};
  end
  if ~cg_params.b_keep_b
      sz = cg_params.NYs(i);  
      bs{i} = x(x_i:x_i+sz-1);
  else
      bs{i} = cg_params.b0{i};
  end
end
