% Generate initial point for SDP
function x = init_sdp(force_marg,bfgs_params,NX,NY,N,px0,py0)

if ~isfield(bfgs_params,'phi0') | isempty(bfgs_params.phi0)
  points = rand(N,NX+NY);
  mn = mean(points,1);
  points = points - repmat(mn,NX+NY,1);
  % Generate their Distance matrix
  x0_mat = my_dist(points,points');
  clear points;
else
  phi0 = bfgs_params.phi0;
  psi0 = bfgs_params.psi0;
  x0 = [phi0;psi0];  
  mn = mean(x0,1);
  x0 = x0 - repmat(mn,NX+NY,1);
  % Generate initial distance matrix
  x0_mat = my_dist(x0,x0');
end

% Translate to PSD n-1 matrix
x0_mat = tauv(x0_mat);
% Add in some noise to the diagonal
x0_mat = x0_mat+eye(NX+NY-1)*1e-4;
% Vectorize
x = smat_to_vec_c(x0_mat)';
clear x0_mat;
% Calculate current normalization factor and make sure starting point is subnormalized
m = vec_to_smat_fast(x,N-1);
dist_m = kappav(m);

if force_marg
  pxy = (px0*py0).*exp(-dist_m(1:NX,NX+1:end));  
else
  pxy = exp(-dist_m(1:NX,NX+1:end));
end

clear dist_m
y = pxy(:);
mxy = max(y);
% The exponent of this is what normalizes the distribution
lse = log(sum(y));

pxy = pxy/exp(lse);

curr_px =  sum(pxy,2);
curr_py =  sum(pxy,1);
max_dev = max(log(curr_px./px0));
max_dev = max(max_dev,max(log(curr_py./py0)));
if max_dev<0
  max_dev = 0;
end
% Add a=log(Z)

if force_marg
  x = [x; lse+max_dev+1e-3];
else
  x = [x; lse+1e-3];
end