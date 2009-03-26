function m = le_trace_solve_gd(nxy,NX,NY,tr_wgt)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

nxy = nxy/sum(nxy(:));

%F={};
% The F matrices take care of the distance matrix being PSD
N = NX+NY;
NVAR = N*(N+1)/2;
trace_grad =sparse(NVAR,1);
x0 = trace_grad;

F = sparse(N*N,NVAR);

uind = triu(ones(N,N));
uind = find(uind');
params.uind = uind;

%F = cell(NVAR,1);
trind = 0;

for i=1:NX+NY
    for j=i:NX+NY
        trind = trind+1;
        mat = sparse(NX+NY,NX+NY);
        mat(i,j) = 1;
        mat(j,i) = 1;
        %        F(:,trind) = mat(:);
        F(:,trind) = mat(:);
        % Generate trace gradient
        if (j==i)
            trace_grad(trind)=tr_wgt;
            x0(trind) = 1;
        end
    end
end

t = 1;
MAXITERS = 100;
niters = 1;
primal_obj = 0;

% Generate posynomial expressions for f0
NVAR = (NX+NY)*(NX+NY+1)/2;
c = sparse(NVAR,1);
%x0  = uptri_to_x(eye(NX+NY))';
A = sparse(NX*NY,NVAR);
trind = 0;

for xi=1:NX
    for yi=NX+1:NX+NY
        trind = trind+1;
        tovec = sparse(1,NVAR);
        tovec(trans_ind(N,xi,yi)) = 2;
        tovec(trans_ind(N,xi,xi)) = -1;
        tovec(trans_ind(N,yi,yi)) = -1;
        
        A(trind,:) = tovec;
        %        c = c+uptri_to_x(mat)*sum(samp(1,:)==xi & samp(2,:)==yi-NX);
        c = c+tovec'*nxy(xi,yi-NX);
    end
end

disp('Done with A');
szs(1) = size(A,1);

if size(szs,1) < size(szs,2),  szs = szs';  end;
n = length(x0);  N = size(A,1);  m = length(szs)-1;
n0 = szs(1);

% E is a matrix s.t. [1'*y0  1'*y1  ... 1'*ym ]' = E*y
E = sparse(m+1,N);
indsf = cumsum(szs)-szs+1;  %  first row of each block of A
indsl= cumsum(szs);  % last row of each block of A      
for i=1:m+1
    E(i,indsf(i):indsl(i)) = ones(1,szs(i));
end;

% INITIALIZE AND CHECK FEASIBILITY OF x0

x = x0;  y = A*x;  f = log(E*exp(y))-dot(c,x);

ntiters = 0;
MAXITERS = 250;  % maximum total number of Newton iterations
NTTOL = 1e-10;    % stop centering when Nt decrement <= NTTOL 
ALPHA = 0.01; 
BETA = 0.5;      % parameters for backtracking line search
MAXLSITERS = 20; % maxiters in line searches backtracking
MU = 2.5;         % mu for update tplus = mu*t;
ABSTOL=1e-7;     % duality gap based stopping criterion
big_iters = 0;
MAXBIGITERS = 15;
CENTER_MINDELTA = 1e-10;
last_f = 0;

params.A = A;
params.c = c;
params.trace_grad = trace_grad;
params.F = F;
params.tr_wgt = tr_wgt ;

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
	x = conj_grad('prev_sdp_grad',params,x0,opts);

    [ft,gt,f] = sdp_grad(x,params);

    if (abs(f-last_f)<CENTER_MINDELTA)
        fprintf(1,'Break Newton at %d\n',center_iters);
        break;
    end
    last_f = f;
    
    t = min(MU*t, (m+1)/ABSTOL);  
    if t>100000
        break;
    end
    big_iters = big_iters+1;
end

m = F*x;   %vec_to_smat(x,F);
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


