function [f,g,basef,newf] = code_gram_sdpgrad(x,params)

global const; % = params.const;
global A; % = params.A;
global px0; % = params.px0;
global py0; % = params.py0;
global xinds; % = params.xinds;
global yinds; % = params.yinds;
global NX;
global NY;
t = params.t;
global c; % = params.c;
global trace_grad;
global tr_wgt;
global logpx0py0;
global replogpx;
global replogpy;
global constlik;

N = NX+NY;

% Create distance matrix in newmat
logZ = x(end);
gram_mat = vec_to_smat_fast(x,N); %F*x;  %
gram_mat = reshape(gram_mat,N,N);

%ev = eig(gram_mat);
[ev,psdness] = chol(gram_mat);
ev = diag(ev);

% Calculate current "partition" function
% Translate to distance matrix
gram_diag = diag(gram_mat);
dist_mat = gram_diag*ones(1,N) + ones(N,1)*gram_diag' - 2*gram_mat;
% Dist mat reduced to [NX,NY]
red_dist_mat2 = dist_mat(NX+1:end,1:NX);

logpxy = logpx0py0' -red_dist_mat2-logZ;
mxy = min(logpxy(:));
% Start out with the log
if max(logpxy(:))>0
  f = Inf;
  g = zeros(size(x));
  return;
else
  log_curr_norm = log(sum(sum(exp(logpxy-mxy))))+mxy;
end

pxy = exp(logpxy);

%log_px_dev  = log(sum(exp(logpxy-mxy-repmat(log(px0'),NY,1)))) + mxy;
%log_py_dev  = log(sum(exp(logpxy-mxy-repmat(log(py0'),1,NX)),2)) + mxy;

log_px_dev  = log(sum(exp(logpxy-mxy-replogpx))) + mxy;
log_py_dev  = log(sum(exp(logpxy-mxy-replogpy),2)) + mxy;


if isinf(log_curr_norm) | isinf(max(abs(log_px_dev))) | ...
isinf(max(abs(log_px_dev))) | log_curr_norm==0 | max(abs(log_px_dev))==0 | max(abs(log_px_dev))==0
     f = Inf;
     g = zeros(size(x));
     return;
end
% Check for violation of constraints
%if (min(ev)<=0) | log_curr_norm>=0  | max(log_px_dev)>=0  | max(log_py_dev)>=0
if (psdness>0) | log_curr_norm>=0  | max(log_px_dev)>=0  | max(log_py_dev)>=0
  f = Inf;
%  fprintf('PSD=%g LogNorm=%g LogPx=%g LogPy=%g\n',psdness>0,log_curr_norm, max(log_px_dev),max(log_py_dev));
  g = zeros(size(x));
  return;
end

logdet = 2*sum(log(ev));

basef = dot(x,c) + tr_wgt*trace(gram_mat);
% Add in all the barriers
%fprintf('Sumall %g\n',sum(gram_mat(:)));

%f = t*basef - sum(log(ev))-sum(-log_px_dev)-sum(-log_py_dev)-log(-log_curr_norm);    %det(gram_mat));     
f = (constlik+basef)+ 1/t*( - logdet -log(-log_curr_norm)-sum(log(-log_px_dev))-sum(log(-log_py_dev)));    %det(gram_mat));     


if nargout==1
  return;
end

bar_grad = det_barrier_grad_gd(gram_mat);

Zgrad2 = sparse(size(dist_mat,1),size(dist_mat,2));
Zgrad2(NX+1:end,1:NX) = 2*pxy;
xsum = 0.5*sum(Zgrad2,2);
ysum = 0.5*sum(Zgrad2,1);
Zgrad2 = Zgrad2 - diag(xsum);
Zgrad2 = Zgrad2 - diag(ysum);
% Add derivative w.r.t logZ
Zgrad2 = smat_to_vec_c(full(Zgrad2));
Zgrad2 = Zgrad2 / exp(log_curr_norm);
Zgrad2 = [Zgrad2 -1];
Zgrad2 = Zgrad2 / log_curr_norm;

tmp_grad = Zgrad2;
xygrad = Zgrad2;

if 0
[i,j] = ndgrid(1:NX,1:NY);
i = i(:); j=j(:);
A = get_gram_op_inds(N,i,NX+j);
A = -A;
A(:,end+1) = -1;


distm = A*x+logpx0py0(:);
y = A*x+logpx0py0(:);
distm = reshape(distm,[NX,NY]);
distm = distm';

pxy2 = exp(distm);
mxy = max(y);
% Start out with the log
norm_expy = y-mxy-log(sum(exp(y-mxy)));
norm_expy = exp(norm_expy);
lse = log(sum(exp(y-mxy)))+mxy;      
gradf = norm_expy;
%      gradf = exp(y)./sum(exp(y));
big_lse = lse;


gradphi =  - 1/lse*A' * gradf; 

for xi=1:NX
  xinds(xi,:) = find(i==xi)';
end
for yi=1:NY
  yinds(yi,:) = find(j==yi)';
end

bigvecx = sparse(size(A,1),1);

% Generate gradient and Hessians for px0 constraints
for xi=1:NX
  yx = y(xinds(xi,:))-log(px0(xi));
  mxy = max(yx);
  % Start out with the log
  lse = log(sum(exp(yx-mxy)))+mxy;      
  norm_expy = yx-lse;
  norm_expy = exp(norm_expy);
  gradf = norm_expy;
  px_lse(xi) = lse;
  bigvecx(xinds(xi,:)) = 1/lse*gradf;
  tmpx{xi} = 1/lse*gradf;
end

gradphix =   - A'*bigvecx;

bigvecy = sparse(size(A,1),1);

for yi=1:NY
  yy = y(yinds(yi,:))-log(py0(yi));
  mxy = max(yy);
  lse = log(sum(exp(yy-mxy)))+mxy;      
  % Start out with the log
  norm_expy = yy-lse; 
  norm_expy = exp(norm_expy);
  gradf = norm_expy;
  py_lse(yi)= lse;
  bigvecy(yinds(yi,:)) = 1/lse*gradf; 
end

gradphiy =  - A'*bigvecy;

end


xsum = sum(pxy,1);
ysum = sum(pxy,2);

pxy_cond = -pxy*spdiags((1./(log_px_dev.*xsum))',0,NX,NX);
newdiag = [-1./log_px_dev';sum(pxy_cond,2) ];
xgrad = sparse(size(dist_mat,1),size(dist_mat,2));  
xgrad(NX+1:end,1:NX) = -2*pxy_cond;
xgrad = xgrad + spdiags(newdiag,0,N,N);

pxy_cond = -spdiags((1./(log_py_dev.*ysum)),0,NY,NY)*pxy;
newdiag = [sum(pxy_cond,1) -1./log_py_dev'];
ygrad = sparse(size(dist_mat,1),size(dist_mat,2));  
ygrad(NX+1:end,1:NX) = -2*pxy_cond;
ygrad = ygrad + spdiags(newdiag',0,N,N);



bar_grad(end+1) = 0;

tmp_grad = xygrad + [smat_to_vec_c(full(xgrad)) -sum(1./log_px_dev)] ...
    + [smat_to_vec_c(full(ygrad)) -sum(1./log_py_dev)];

g =(c'+trace_grad)+ 1/t*(bar_grad -tmp_grad') ; 

