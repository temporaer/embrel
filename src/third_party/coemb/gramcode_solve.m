function [m,train_lik] = gramcode_solve(nxy,NX0,NY0,tr_wgt0)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

global const; 
global px0; 
global py0; 
global NX; 
global NY; 
global c; 
global trace_grad;
global tr_wgt;
global logpx0py0;
global replogpx;
global replogpy;
global constlik;

NX = size(nxy,1);
NY = size(nxy,2);
N = NX+NY;

% Normalize to the number of eigenvalues
tr_wgt = tr_wgt0; %/N;

nxy = nxy/sum(nxy(:));

px0 = sum(nxy,2);
py0 = sum(nxy,1);
replogpx  = repmat(log(px0'),NY,1);
replogpy  = repmat(log(py0'),1,NX);


logpx0py0 = log(px0*py0);

%F={};
% The F matrices take care of the distance matrix being PSD

NM1 = N-1;
NVAR = (N-1)*N/2;
trace_grad =sparse(NVAR,1);

% Initialize x0 to NX+NY points radnomly spaced around the origin
rand('seed',0);

% Generate the vector c, for the linear target function
mat = sparse(N,N);
mat(1:NX,NX+1:end) = nxy;
px = sum(mat,2);
py = sum(mat,1);
mat = mat + diag(px);
mat = mat + diag(py);
mat(1:NX,NX+1:end) = -2*mat(1:NX,NX+1:end);
%mat(NX+1:end,1:NX) = mat(1:NX,NX+1:end)';
c = smat_to_vec_c(full(mat'));
% Add an element for +log(Z)
c(end+1) = 1;

constlik = -dot(logpx0py0(:),nxy(:));

diag_inds = smat_to_vec_c(eye(N))';
trace_grad = diag_inds*tr_wgt;
trace_grad(end+1) = 0;

% Initialize x0 to NX+NY points radnomly spaced around the origin
rand('seed',0);
points = rand(N,NX+NY)*1e-3;
%mn = mean(points,1);
% Generate their Gram matrix
x0_mat = points*points';
%x0_mat = eye(N)*1e-3;
% Vectorize
x0 = smat_to_vec_c(x0_mat)';

% Calculate current normalization factor and make sure starting point is subnormalized
dist_m = g_to_d(x0_mat);
pxy = (px0*py0).*exp(-dist_m(1:NX,NX+1:end));
y = pxy(:);
mxy = max(y);
% Start out with the log
lse = log(sum(exp(y-mxy)))+mxy      
% Add a=log(Z)


x0 = [x0; lse];

t = 1;
MAXITERS = 1000;
niters = 1;
primal_obj = 0;

% Generate posynomial expressions for f0

%E = ones(1,size(A,1));

% INITIALIZE AND CHECK FEASIBILITY OF x0


x = x0;  f =dot(c,x);

last_train_lik = Inf;
ntiters = 0;
MAXITERS = 250;  % maximum total number of Newton iterations
NTTOL = 1e-3;    % stop centering when Nt decrement <= NTTOL 
ALPHA = 0.01; 
BETA = 0.5;      % parameters for backtracking line search
MAXLSITERS = 20; % maxiters in line searches backtracking
MU = 10;         % mu for update tplus = mu*t;
ABSTOL=1e-7;     % duality gap based stopping criterion
big_iters = 0;
MAXBIGITERS = 15;
CENTER_MINDELTA = 1e-5;
DIFFTOL = 1e-3;
last_f = 0;
train_lik = 0;

% Correct this
m = NVAR;

params.t = t;

while (big_iters <= MAXBIGITERS)
  
  loglik = 0;
  fprintf(1,'t=%g  N=%4d f=%g. Lik=%g\n',t,big_iters,full(f),full(train_lik));   
  
  x0 = x;
  params.t = t;
%code_gram_sdpgrad'  
  x = bfgs_mine(x0,'code_gram_sdpgrad',1e-6,1e5,0,params);
  train_lik = dot(c,x);
  
  t = MU*t;
  rt =   abs((last_train_lik-train_lik)/last_train_lik);
  fprintf('Rat=%g\n',rt);
  if rt<DIFFTOL
     break;
  end
  last_train_lik = train_lik;
  
  big_iters = big_iters+1;
end

% This is (n-1) SDP
m = vec_to_smat_fast(x(1:end-1),N);
% Transform it into distance matrix
k = NX+NY;
m = reshape(m,k,k);
dist_m = g_to_d(m);
% Translate to centered Gram
V = eye(k)-1/k*ones(k,1)*ones(1,k);
cent_m = -0.5*V*dist_m*V;

sq = sqrtm(full(cent_m));
sq  = sq-repmat(mean(sq),size(sq,1),1);
m  = sq*sq'; 

pxy = (px0*py0).*exp(-dist_m(1:NX,NX+1:end)-x(end));
y = pxy(:);
mxy = max(y);
% Start out with the log
lse = log(sum(y));      

fprintf('NormTo=%g\n',exp(lse));
fprintf('Mypy: %s\n',num2str(sum(pxy,1)));
fprintf('py0: %s\n',num2str(py0));
fprintf('Mypx: %s\n',num2str(sum(pxy,2)'));
fprintf('px0: %s\n',num2str(px0'));


function v = trans_ind(N,xi,yi)

v = sum(N:-1:N-xi+2)+(yi-xi+1);


