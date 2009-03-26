function [m,train_lik] = code_sdp_solve_grad(nxy,NX0,NY0,tr_wgt0)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

global const; 
global A; 
global px0; 
global py0; 
global xinds; 
global yinds; 
global NX; 
global NY; 
global c; 
global trace_grad;
global tr_wgt;

tr_wgt = tr_wgt0;

nxy = nxy/sum(nxy(:));

px0 = sum(nxy,2);
py0 = sum(nxy,1);

%F={};
% The F matrices take care of the distance matrix being PSD
NX = size(nxy,1);
NY = size(nxy,2);
N = NX+NY;
NM1 = N-1;
NVAR = (N-1)*N/2;
trace_grad =sparse(NVAR,1);

% Initialize x0 to NX+NY points radnomly spaced around the origin
rand('seed',0);

% Find how all the dxy are generated from p
[i,j] = ndgrid(1:NX,1:NY);
i = i(:); j=j(:);
A = get_kappa_op_inds(N,i,NX+j,ones(1,length(i)));

[nxyinds,tmp,v] = find(nxy(:));
% Translate to ot
c = v'*A(nxyinds,:);
c = c';
c(end+1) = 1;


A = -A;
%[xi,yi,v] = find(nxy);
% Get gradient of the \sum n(x,y)d(x,y) element
%[c] = get_kappa_op_inds(N,xi,NX+yi,v);
%c = c';
% This element is for +log(Z)
%c(end+1) = 1;
% This is for the condition \sum_{x,y} exp(-d(xy)-a) <=1
A(:,end+1) = -1;

xinds = zeros(NX,NY);
for xi=1:NX
  xinds(xi,:) = find(i==xi)';
end
yinds = zeros(NY,NX);
for yi=1:NY
  yinds(yi,:) = find(j==yi)';
end

for k=1:length(i)
  const(k,1) = log(px0(i(k)))+log(py0(j(k)));
end


% Get gradient of the trace element
%trace_grad = get_kappa_op_inds(N,1:N,1:N,ones(1,N));
DIFFTOL = 1e-3; 


F = sparse((N-1)^2,NVAR);
F = cell(NVAR,1);
trind = 0;

lineq_vec = zeros(1,NVAR);

for i=1:N-1
  for j=i:N-1
    trind = trind+1;
    mat = sparse(N-1,N-1);
    mat(i,j) = 1;
    mat(j,i) = 1;
    %        F(:,trind) = mat(:);
    F{trind} = mat;
    if (j==i)
      trace_grad(trind)=tr_wgt;
    end
  end
end

trace_grad(end+1,1)=0;

% Initialize x0 to NX+NY points radnomly spaced around the origin
rand('seed',0);
points = rand(N,NX+NY);
mn = mean(points,1);
points = points - repmat(mn,NX+NY,1);
% Generate their Distance matrix
x0_mat = my_dist(points,points');
% Translate to PSD n-1 matrix
x0_mat = tauv(x0_mat);
% Vectorize
x0 = smat_to_vec_c(x0_mat)';

% Calculate current normalization factor and make sure starting point is subnormalized
m = vec_to_smat(x0,F);
dist_m = kappav(m);
pxy = (px0*py0).*exp(-dist_m(1:NX,NX+1:end));
y = pxy(:);
mxy = max(y);
% Start out with the log
lse = log(sum(exp(y-mxy)))+mxy      
% Add a=log(Z)


x0 = [x0; 20];

t = 1;
niters = 1;
primal_obj = 0;

% Generate posynomial expressions for f0


% INITIALIZE AND CHECK FEASIBILITY OF x0


x = x0;  f =dot(c,x);


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
A_lineq = lineq_vec;
n_lineq = size(A_lineq,1);

% This is the order of constraints we have, which determines the duality gap 
% for a given t value

m = 2*(NX+NY);

params.t = t;

while (big_iters <= MAXBIGITERS)
  
  loglik = 0;
  print_to_log('t=%g  N=%d f=%g. Lik=%g\n',t,big_iters,full(f),full(train_lik));   
  
  x0 = x;
  params.t = t;
  x = bfgs_mine(x0,'code_sdp_grad',1e-6,1e5,0,params);

  train_lik = dot(c,x);
  newmat = vec_to_smat(x(1:end-1),F);
  f = train_lik+tr_wgt*trace(newmat);
  
%  t = min(MU*t, (m+1)/ABSTOL);  
  t = MU*t;
%  if t>(m+1)/ABSTOL
%    break;
%  end
  rt =   abs((last_train_lik-train_lik)/last_train_lik);

  if rt<DIFFTOL
     break;
  end
  last_train_lik = train_lik;
  big_iters = big_iters+1;
end

% This is (n-1) SDP
m = vec_to_smat(x(1:end-1),F);
% Transform it into distance matrix
dist_m = kappav(m);

% Translate to centered Gram
k = NX+NY;
V = eye(k)-1/k*ones(k,1)*ones(1,k);
cent_m = -0.5*V*dist_m*V;

sq = sqrtm(full(cent_m));
sq  = sq-repmat(mean(sq),size(sq,1),1);
m  = sq*sq'; %m-repmat(mean(m),size(m,1),1);

%for k=0:length(x)-1
%    v = [v diag(x,k)];
%end

pxy = (px0*py0).*exp(-dist_m(1:NX,NX+1:end)-x(end));
y = pxy(:);
mxy = max(y);
% Start out with the log
lse = log(sum(y));      

print_to_log('NormTo=%g\n',exp(lse));
print_to_log('Mypy: %s\n',num2str(sum(pxy,1)));
print_to_log('py0: %s\n',num2str(py0));
print_to_log('Mypx: %s\n',num2str(sum(pxy,2)'));
print_to_log('py0: %s\n',num2str(px0'));


function v = trans_ind(N,xi,yi)

v = sum(N:-1:N-xi+2)+(yi-xi+1);


