function [m,train_lik] = cond_sdp_code(nxy,NX0,NY0,tr_wgt0,bfgs_params)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

global px0; 
global py0; 
global xinds; 
global yinds; 
global NX; 
global NY; 
global c; 
global trace_grad;
global tr_wgt;
global pxy0;



tr_wgt = tr_wgt0;

nxy = nxy/sum(nxy(:));

px0 = sum(nxy,2);
py0 = sum(nxy,1);
pxy0 = nxy;

%F={};
% The F matrices take care of the distance matrix being PSD
NX = size(nxy,1);
NY = size(nxy,2);
N = NX+NY;
NM1 = N-1;
NVAR = (N-1)*N/2;
trace_grad =sparse(NVAR,1);

c = make_c(pxy0,NX,NY,N-1);
c = smat_to_vec_c(c')';
c(end+1:end+NX) = px0;

[i,j] = ndgrid(1:NX,1:NY);
i = i(:); j=j(:);

xinds = zeros(NX,NY);
for xi=1:NX
  xinds(xi,:) = find(i==xi)';
end
yinds = zeros(NY,NX);
for yi=1:NY
  yinds(yi,:) = find(j==yi)';
end
clear i j 

% Get gradient of the trace element
%trace_grad = get_kappa_op_inds(N,1:N,1:N,ones(1,N));
DIFFTOL = 1e-6; 

z = eye(N-1);
trace_grad = smat_to_vec(z)*tr_wgt;
trace_grad(1,end+1:end+NX)=0;
trace_grad = trace_grad';
clear z;

% Initialize x0 to NX+NY points radnomly spaced around the origin
rand('seed',0);
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

pxy = exp(-dist_m(1:NX,NX+1:end)+repmat(log(py0),NX,1));

xlse = log(sum(pxy,2));
clear dist_m

xlse(xlse<0) = 0;

% Add a=log(Z)

x = [x; xlse+1e-4];

t = bfgs_params.t0;
niters = 1;
primal_obj = 0;

% Generate posynomial expressions for f0
f =dot(c,x);


ntiters = 0;
MU = 10; %2.5;         % mu for update tplus = mu*t;
ABSTOL=1e-7;     % duality gap based stopping criterion
big_iters = 0;
MAXBIGITERS = 15;
last_f = 0;
train_lik = 0;
last_train_lik = Inf;
% Generate equality constraint which says matrix elements should
% sum to zero
%A_lineq = lineq_vec;
%n_lineq = size(A_lineq,1);

% This is the order of constraints we have, which determines the duality gap 
% for a given t value

m = 2*(NX+NY);

params.t = t;

fun = 'cond_sdp_grad';
  
% INITIALIZE AND CHECK FEASIBILITY OF x0
f0 = feval(fun,x,params);
if isinf(f0)
  disp('Infeasible start');
  return;
else
  disp('Feasible start');
end

if isfield(bfgs_params,'tmp_fname') & xistf(bfgs_params.tmp_fname)
   load(bfgs_params.tmp_fname);
   print_to_log('Commencing at t=%g\n',t);
end

while (big_iters <= MAXBIGITERS)
  
  loglik = 0;
  print_to_log('t=%g  N=%4d f=%g. Lik=%g\n',t,big_iters,full(f),full(train_lik));   
  
  if isfield(bfgs_params,'tmp_fname') 
	save(bfgs_params.tmp_fname); 
  end

  params.t = t;
  tol = 1e-6; %bfgs_params.bfgs_tol;
  bfgs_niter = 250; %bfgs_params.bfgs_niter;
  bfgs_nsmax = bfgs_params.bfgs_nsmax;

  x = bfgs_mine(x,fun,tol,bfgs_niter,0,params,bfgs_nsmax);    
  
  train_lik = dot(c,x);
  newmat = vec_to_smat_fast(x(1:end-NX),N-1);
  f = train_lik+tr_wgt*trace(newmat);

  newmat=[];
  
  t = MU*t;
  rt =   abs((last_train_lik-train_lik)/last_train_lik);
  if rt<DIFFTOL
     break;
  end
  last_train_lik = train_lik;
  big_iters = big_iters+1;
end

% This is (n-1) SDP
m = vec_to_smat_fast(x(1:end-NX),N-1);
% Transform it into distance matrix
dist_m = kappav(m);

% Translate to centered Gram
k = NX+NY;
V = eye(k)-1/k*ones(k,1)*ones(1,k);
cent_m = -0.5*V*dist_m*V;

sq = sqrtm(full(cent_m));
sq  = sq-repmat(mean(sq),size(sq,1),1);
m  = sq*sq'; %m-repmat(mean(m),size(m,1),1);

zxs = x(end-NX+1:end); 
pxy = exp(-dist_m(1:NX,NX+1:end)-repmat(zxs,1,NY)+repmat(log(py0),NX,1));
xlse = full(sum(pxy,2)); 

print_to_log('NormTo=%g-%g\n',min(xlse),max(xlse));


function v = trans_ind(N,xi,yi)

v = sum(N:-1:N-xi+2)+(yi-xi+1);


