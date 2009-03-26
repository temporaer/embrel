function m = le_hess(nxy,NX,NY,tr_wgt)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

nxy = nxy/sum(nxy(:));

%F={};
% The F matrices take care of the distance matrix being PSD
N = NX+NY;
NVAR = N*(N+1)/2;
trace_grad =sparse(NVAR,1);
x0 = trace_grad;

trind = 0;

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
    end
end

mat = sparse(N,N);
mat(1:NX,NX+1:end) = nxy;
px = sum(mat,2);
py = sum(mat,1);
mat = mat + diag(px);
mat = mat + diag(py);
mat(1:NX,NX+1:end) = -2*mat(1:NX,NX+1:end);
mat(NX+1:end,1:NX) = mat(1:NX,NX+1:end)';
c = smat_to_vec_c(full(mat))';

E = ones(1,size(A,1));

% INITIALIZE AND CHECK FEASIBILITY OF x0
x0 = smat_to_vec_c(eye(N))';
x = x0;  y = A*x;  f = log(sum(exp(y)))-dot(c,x);

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
CENTER_MINDELTA = 1e-5;
last_f = 0;

while (big_iters <= MAXBIGITERS)

   % Next lines for explicitly calculating target function   
   loglik = 0;
   fprintf(1,'t=%g  N=%4d f=%g. RealLog=%g\n',t,big_iters,full(f),loglik);   
   
   center_iters = 0;
   while (center_iters <= MAXITERS)

     % FORM GRADIENT AND HESSIAN
     
     gram_mat = vec_to_smat_fast(x,N);     
     dist_mat = gram_diag*ones(1,N) + ones(N,1)*gram_diag' - 2*gram_mat;
     red_dist_mat = dist_mat(NX+1:end,1:NX);
     
     [bar_grad,bar_hess] = det_barrier_hess(gram_mat);
     logZ = log(sum(exp(-red_dist_mat2(:)+red_dist_mat2(1,1))))-red_dist_mat(1,1);         
     
     % Generate function for the exponent
     %      gradf = exp(y)./(E'*(E*exp(y)));
 
if 0      
     gradf = exp(y)./sum(exp(y));
   
     gradphi = A' * (gradf*t)-t*c +bar_grad + t*trace_grad; 
     hessphi =  - ... % add the rank one hessf terms
	 (A' * (sparse(diag(exp(y))) * E' * ...
	 t./(E*exp(y)).^2)) * ...
	 ((E * sparse(diag(exp(y)))) * A)  ... 
	 + ... % add the diagonal hessf terms
	 (A'*sparse(diag(gradf.*(E'* t)))) * A + bar_hess; % + centering_hess;
   
end

     off_diag = sparse(size(dist_mat,1),size(dist_mat,2));
     off_diag(NX+1:end,1:NX) = exp(-red_dist_mat-logZ);
     xsum = sum(off_diag,2);
     ysum = sum(off_diag,1);
     off_diag(NX+1:end,1:NX) = 2*off_diag(NX+1:end,1:NX);
     
     on_diag = sparse(size(dist_mat,1),size(dist_mat,2));     
     on_diag = on_diag - diag(xsum);
     on_diag  = on_diag - diag(ysum);
     
     Zgrad = smat_to_vec_c(full(off_diag+on_diag));
     
     gradphi = t*(Zgrad'+c'+trace_grad)+ bar_grad  ; 
     
     
     
     if (rcond(full(hessphi))<eps)
       fprintf(1,'Ill conditioned in Newton\n');
       break;
     end
     dx = -hessphi\gradphi;
     %      dx = [0;dx];       % Generate zero movement on the first element
     center_iters = center_iters + 1;
     %      fprime = gradphi'*dx;
     dy = A*dx;
     %      if (-fprime/2 < NTTOL), break; end; 
     
     
     % BACKTRACKING LINE SEARCH
     

     phi = t*f - log(det(gram_mat)); %-log(sum(xmat)); %sum(log(-f(2:m+1)));
     s = 1e-1;
     for lsiters = 1:MAXLSITERS
       newmat = vec_to_smat_fast(x+s*dx,N);
       
       newf = log(E*exp(y+s*dy))-dot(c,(x+s*dx)) + tr_wgt*trace(newmat);
       
       if (min(eig(newmat))> 0) % & sum(newmat)>0
	 newphi = t*newf - log(det(newmat)); 
	 %            if (newphi < phi + ALPHA*s*fprime)
	 if (newphi < phi )
	   break;
	 end
       end
       s = BETA*s;
     end
     
     x = x+s*dx;
     y = A*x;
     newmat = vec_to_smat_fast(x+s*dx,N);
     f = log(E*exp(y))-dot(c,x)+tr_wgt*trace(newmat);
     %      fprintf('The trace is %g\nEigs:',trace(newmat));
     %      eig(newmat)
     if (abs(f-last_f)<CENTER_MINDELTA)
       fprintf(1,'Break Newton at %d\n',center_iters);
       break;
     end
     last_f = f;
     
     %     fprintf('In Newton %g. Iter=%d\n',newphi,lsiters);
   end
   
   
   % DUAL VARIABLES ON CENTRAL PATH
   %
   % nu0 = (1 - g0'*dy0)*g0 + diag(dy0)*g0 
   % nui = (-1/t*fi)*( (1-(1+1/fi)*gi'*dyi)*gi + diag(dyi)*gi )
   % lambdai = (-1/tfi)*(1-gi'*dyi/fi) 
   %
   % these expressions are accurate even with imperfect centering
   %


%   primalobj = f(1);
%   gap = primalobj-dualobj; 

   t = min(MU*t, (m+1)/ABSTOL);  
   if t>100000
       break;
   end
   big_iters = big_iters+1;
end

m = vec_to_smat_fast(x,N);
dist_m = g_to_d(m);
k = NX+NY;
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


if 0
       off_diag = sparse(size(dist_mat,1),size(dist_mat,2));
     off_diag(NX+1:end,1:NX) = exp(-red_dist_mat-logZ);
     xsum = sum(off_diag,2);
     ysum = sum(off_diag,1);
     off_diag(NX+1:end,1:NX) = 2*off_diag(NX+1:end,1:NX);
     
     on_diag = sparse(size(dist_mat,1),size(dist_mat,2));     
     on_diag = on_diag - diag(xsum);
     on_diag  = on_diag - diag(ysum);

     Zgrad = smat_to_vec_c(full(off_diag+on_diag));

     objective_grad = t*(Zgrad'+c'+trace_grad) ; 
end