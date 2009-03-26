function [f,g] = mway_logemb_grad(x,params)

[phi,psis,as,bs] = read_mway_params(x,params);

f = 0;
phigrad =0;
ndist = length(params.pxy0s);

for i=1:ndist
  [val,gphi,gpsi{i},ga{i},gb{i}] = ...
      mrgle_grad_fast_prm(phi,psis{i},as{i},bs{i},sparse(params.pxy0s{i}),params.px0s{i},params.py0s{i},params.b_cond);

  invpx = spdiags(1./params.px0s{i},0,length(params.px0s{i}),length(params.px0s{i}));
  invpy = spdiags(1./params.py0s{i}',0,length(params.py0s{i}),length(params.py0s{i}));
  
%  invpx = diag(1./params.px0s{i});
%  invpy = diag(1./params.py0s{i});

  ga{i} = ga{i}*params.w_nxy(i);
  
  gb{i} = gb{i}*params.w_nxy(i);

  phigrad = phigrad+params.w_nxy(i)*gphi;
  f = f + params.w_nxy(i)*val;
end

if isfield(params,'pxx0') & ~isempty(params.pxx0)
   pxx_x0 = sum(params.pxx0,2);
  [val,gphi,dummy,tmpga] = ...
      mrgle_grad_fast_prm(phi,phi,as{1},as{1},sparse(params.pxx0),pxx_x0,pxx_x0,params.b_cond);
   ga{1} = ga{1} + 2*tmpga*params.w_pxx0;  
   f = f + params.w_pxx0*val;
   phigrad = phigrad + params.w_pxx0*gphi;
end


if isfield(params,'pyy0') & ~isempty(params.pyy0)
   pyy_y0 = sum(params.pyy0,1);
  [val,tmp_gpsi,dummy,tmpgb] = ...
      mrgle_grad_fast_prm(psis{1},psis{1},bs{1},bs{1},sparse(params.pyy0),pyy_y0,pyy_y0,params.b_cond);
   gb{1} = gb{1} + 2*tmpgb*params.w_pyy0;  
   f = f + params.w_pyy0*val;
   gpsi{1} = gpsi{1} + params.w_pyy0*tmp_gpsi;
end

% Construct overall gradient
g  = [phigrad(:)];
if params.b_keep_phi
  g = [];
else
  g = [phigrad(:)];
end

for i=1:length(params.pxy0s)
  if ~params.b_keep_psi
    g = [g;gpsi{i}(:)];
  end
  
  if ~params.b_keep_a
    g = [g;ga{i}(:)]; %ga{i} = zeros(size(ga{i}));
  end
  
  if ~params.b_keep_b
    g = [g;gb{i}(:)]; %gb{i} = zeros(size(gb{i}));
  end

%  g = [g;gpsi{i}(:);ga{i}(:);gb{i}(:)];
end

