function m = le_trace_solve_fixorig(nxy,NX,NY,tr_wgt)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

nxy = nxy/sum(nxy(:))

F={}; F2={}; 
% The F matrices take care of the distance matrix being PSD
N = NX+NY-1; % Total number of points
NVAR = N*(N+1)/2;
trace_grad =zeros(NVAR,1);
stable_inds=[];   %Indices where the y_NY element appears. They are identically zero

for i=1:N
    for j=i:N
        mat = zeros(N,N);
        mat(i,j) = 1;
        mat(j,i) = 1;
        F{end+1} = mat;
        mat = zeros(N+1,N+1);
        mat(i,j) = 1;
        mat(j,i) = 1;
        F2{end+1} = mat;
        % Generate trace gradient
        if (j==i)
            trace_grad(length(F))=tr_wgt;
        end
    end
end

t = 1;
MAXITERS = 100;
niters = 1;
primal_obj = 0;

% Generate posynomial expressions for f0
c = zeros(1,NVAR)
x0  = uptri_to_x(eye(N))';
A=[]; b= [];
for xi=1:NX
    for yi=NX+1:NX+NY
        mat = zeros(N,N);
        if yi~=NX+NY
            mat(xi,yi) = 2;
            mat(yi,yi)=-1;
        end
        mat(xi,xi)=-1;        
        A(end+1,:) = uptri_to_x(mat);
        b(end+1) = 0;
%        c = c+uptri_to_x(mat)*sum(samp(1,:)==xi & samp(2,:)==yi-NX);
        c = c+uptri_to_x(mat)*nxy(xi,yi-NX);
    end
end
b = b';
%c = c'/length(samp);
c = c'/sum(nxy(:));

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

x = x0; 
y = A*x+b;  f = log(E*exp(y))-dot(c,x);

ntiters = 0;
MAXITERS = 40;  % maximum total number of Newton iterations
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
   curr_dmat = g_to_d(vec_to_smat(x,F2));
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
   
   fprintf(1,'t=%g  N=%4d f=%g. RealLog=%g\n',t,big_iters,f,loglik);   
   
   center_iters = 0;
   while (center_iters <= MAXITERS)

      % FORM GRADIENT AND HESSIAN

      [bar_grad,bar_hess] = det_barrier_grad(x,F);

      % Generate function for the exponent
%      gradf = exp(y)./(E'*(E*exp(y)));
      gradf = exp(y)./sum(exp(y));
      
      xmat = vec_to_smat(x,F);     
      
      gradphi = A' * (gradf*t)-t*c +bar_grad + trace_grad; 
      hessphi =  - ... % add the rank one hessf terms
                (A' * (sparse(diag(exp(y))) * E' * ...
                t./(E*exp(y)).^2)) * ...
                ((E * sparse(diag(exp(y)))) * A)  ... 
              + ... % add the diagonal hessf terms
                (A'*sparse(diag(gradf.*(E'* t)))) * A + bar_hess; % + centering_hess;

      if (rcond(hessphi)<eps)
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
         newmat = vec_to_smat(x+s*dx,F);
         
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
      y = A*x+b;
      newmat = vec_to_smat(x+s*dx,F);
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

  if (abs(f-last_f)>=CENTER_MINDELTA)
      fprintf(1,'Need more  Newtons at %d\n',center_iters);
      break;
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

m = vec_to_smat(x,F2);
dist_m = g_to_d(m);
k = NX+NY;
V = eye(k)-1/k*ones(k,1)*ones(1,k);
cent_m = -0.5*V*dist_m*V;

sq = sqrtm(m);
sq  = sq-repmat(mean(sq),size(sq,1),1);
m  = sq*sq'; %m-repmat(mean(m),size(m,1),1);

eig(vec_to_smat(x,F))
%for k=0:length(x)-1
%    v = [v diag(x,k)];
%end



