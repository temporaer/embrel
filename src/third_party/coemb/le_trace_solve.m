function [m,train_lik] = le_trace_solve(nxy,NX,NY,tr_wgt)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

nxy = nxy/sum(nxy(:));

%F={};
% The F matrices take care of the distance matrix being PSD
N = NX+NY;
NVAR = N*(N+1)/2;
trace_grad =sparse(NVAR,1);
x0 = trace_grad;

F = sparse(N*N,NVAR);
F = cell(NVAR,1);
trind = 0;

for i=1:NX+NY
    for j=i:NX+NY
        trind = trind+1;
        mat = sparse(NX+NY,NX+NY);
        mat(i,j) = 1;
        mat(j,i) = 1;
%        F(:,trind) = mat(:);
        F{trind} = mat;
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
CENTER_MINDELTA = 1e-5;
last_f = 0;
train_lik = 0;

while (big_iters <= MAXBIGITERS)

   % Next lines for explicitly calculating target function   
   if 0 
   curr_dmat = g_to_d(vec_to_smat(x,F));
   curr_dmat = curr_dmat(1:NX,NX+1:end);
   pxy = exp(-curr_dmat);
   Z = sum(pxy(:));
%   logpxy = -curr_dmat-log(Z);  %pxy/sum(pxy(:));
%   loglik = 0;
    logpxy = -curr_dmat-log(Z);
%   for i=1:length(samp)
%       loglik = loglik + logpxy(samp(1,i),NX+samp(2,i));
%   end
   loglik = sum(sum(nxy.*logpxy));
   loglik = loglik/length(sum(sum(nxy)));
else
    loglik = 0;
end
   fprintf(1,'t=%g  N=%4d f=%g. Lik=%g\n',t,big_iters,full(f),full(train_lik));   
   
   center_iters = 0;
   while (center_iters <= MAXITERS)

      % FORM GRADIENT AND HESSIAN

      [bar_grad,bar_hess] = det_barrier_grad(x,F);

      % Generate function for the exponent
%      gradf = exp(y)./(E'*(E*exp(y)));
      gradf = exp(y)./sum(exp(y));
      
      xmat = vec_to_smat(x,F);     
      
      gradphi = A' * (gradf*t)-t*c +bar_grad + t*trace_grad; 
      hessphi =  - ... % add the rank one hessf terms
                (A' * (sparse(diag(exp(y))) * E' * ...
                t./(E*exp(y)).^2)) * ...
                ((E * sparse(diag(exp(y)))) * A)  ... 
              + ... % add the diagonal hessf terms
                (A'*sparse(diag(gradf.*(E'* t)))) * A + bar_hess; % + centering_hess;

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


      phi = t*f - log(det(xmat)); %-log(sum(xmat)); %sum(log(-f(2:m+1)));
      s = 1e-1;
      for lsiters = 1:MAXLSITERS
         newmat = full(vec_to_smat(x+s*dx,F));
         
         newf = log(E*exp(y+s*dy))-dot(c,(x+s*dx)) + tr_wgt*trace(newmat);

         if (min(eig(newmat))> 0) % & sum(newmat)>0
            newphi = t*newf - log(det(newmat)); 
            if (newphi < phi + ALPHA*s*fprime)
%            if (newphi < phi )
              break;
            end
         end
         s = BETA*s;
      end

      x = x+s*dx;
      y = A*x;
      newmat = vec_to_smat(x+s*dx,F);
      train_lik = log(E*exp(y))-dot(c,x);
      f = train_lik+tr_wgt*trace(newmat);
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

m = vec_to_smat(x,F);
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


