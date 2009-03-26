function [m,train_lik] = le_trace_gd_str(nxy,NX,NY,tr_wgt,x0)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

nxy = nxy/sum(nxy(:));

%F={};
% The F matrices take care of the distance matrix being PSD
N = NX+NY;
NVAR = N*(N+1)/2;
trace_grad =sparse(NVAR,1);
c = sparse(1,NVAR);
%x0 = trace_grad;

%F = sparse(N*N,NVAR);

%F = cell(NVAR,1);
trind = 0;

k = [1 N];

uind = 0; %sparse(1,NVAR);

if ~exist('x0')
    x0 = smat_to_vec_c(eye(N))';
else
    x0 = smat_to_vec_c(x0)';
end

%x0 = smat_to_vec_c(diag(rand(1,N)))';
trace_grad = x0*tr_wgt;

params.uind = uind;

mat = sparse(N,N);
mat(1:NX,NX+1:end) = nxy;
px = sum(mat,2);
py = sum(mat,1);
mat = mat + diag(px);
mat = mat + diag(py);
mat(1:NX,NX+1:end) = -2*mat(1:NX,NX+1:end);
%mat(NX+1:end,1:NX) = mat(1:NX,NX+1:end)';
c = smat_to_vec_c(full(mat'));

if 0 
trind = 0;
for xi=1:NX
    for yi=NX+1:NX+NY
        if nxy(xi,yi-NX)==0
	  continue;
        end
        trind = trind+1;
        tovec = sparse(1,NVAR);
        tovec(trans_ind(N,xi,yi)) = -2;
        tovec(trans_ind(N,xi,xi)) = 1;
        tovec(trans_ind(N,yi,yi)) = 1;
        c = c+tovec*nxy(xi,yi-NX);
    end
end
end

t = 1;
MAXITERS = 100;
niters = 1;
primal_obj = 0;

% Generate posynomial expressions for f0
NVAR = (NX+NY)*(NX+NY+1)/2;

%x0  = uptri_to_x(eye(NX+NY))';
% INITIALIZE AND CHECK FEASIBILITY OF x0

ntiters = 0;
MAXITERS = 250;  % maximum total number of Newton iterations
NTTOL = 1e-10;    % stop centering when Nt decrement <= NTTOL 
ALPHA = 0.01; 
BETA = 0.5;      % parameters for backtracking line search
MAXLSITERS = 20; % maxiters in line searches backtracking
MU = 10; %2.5;         % mu for update tplus = mu*t;
ABSTOL=1e-7;     % duality gap based stopping criterion
big_iters = 0;
MAXBIGITERS = 15;
CENTER_MINDELTA = 1e-4;
last_f = 0;
m = 1;

params.trace_grad = trace_grad;
%params.F = F;
params.tr_wgt = tr_wgt ;
params.c = c;
params.nxy = nxy;
params.NX = NX;
params.NY = NY;

last_ftrace = 0;
f = 0;

params.t = 0;
x = x0;
[ft,gt,f] = sdp_grad(x,params);

while (big_iters <= MAXBIGITERS)
    
    % Next lines for explicitly calculating target function   
    loglik = 0;
    fprintf(1,'t=%g  N=%4d f=%g. RealLog=%g\n',t,big_iters,full(f),loglik);   
    
    center_iters = 0;
    % Generate function for the exponent
    center_iters = center_iters + 1;
    
    params.t = t;
    opts = -1*zeros(1,9);
    x0 = x;
    if big_iters>=3
        x = conj_grad('sdp_grad',params,x0,opts);
    else
        x = bfgswopt(x0,'sdp_grad',params,1e-6,1e5);    
    end
    
    [ft,gt,f,ftrace] = sdp_grad(x,params);

    if (abs(f-last_f)<CENTER_MINDELTA)
        fprintf(1,'Break Newton at %d\n',center_iters);
        break;
    end
    last_f = f;
    
    t = min(MU*t, (m+1)/ABSTOL);  
    if t>1000
        break;
    end
    big_iters = big_iters+1;
end

train_lik = full(f);

m = vec_to_smat_fast(x,N);
%F*x;   %vec_to_smat(x,F);
k = NX+NY;
m = reshape(m,k,k);
dist_m = g_to_d(m);
V = eye(k)-1/k*ones(k,1)*ones(1,k);
cent_m = -0.5*V*dist_m*V;

sq = sqrtm(full(m));
sq  = sq-repmat(mean(sq),size(sq,1),1);
m  = sq*sq'; %m-repmat(mean(m),size(m,1),1);

%for k=0:length(x)-1
%    v = [v diag(x,k)];
%end

function v = trans_ind(N,xi,yi)

v = sum(N:-1:N-xi+2)+(yi-xi+1);


