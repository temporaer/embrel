function [m,train_lik] = code_sdp_solve2(nxy,NX,NY,tr_wgt)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

nxy = nxy/sum(nxy(:));

%F={};
% The F matrices take care of the distance matrix being PSD
N = NX+NY;
NM1 = N-1;
NVAR = (N-1)*N/2;
trace_grad =sparse(NVAR,1);

% Initialize x0 to NX+NY points radnomly spaced around the origin
rand('seed',0);
[xi,yi,v] = find(nxy);
% Get gradient of the \sum n(x,y)d(x,y) element
[d_grad] = get_kappa_op_inds(N,xi,NX+yi,v);
d_grad = -d_grad';
c = d_grad;
% Find how all the dxy are generated from p
[i,j] = meshgrid(1:NX,1:NY);
[tmp,A] = get_kappa_op_inds(N,i(:),NX+j(:),ones(1,length(i(:))));
A = -A;
% Get gradient of the trace element
%trace_grad = get_kappa_op_inds(N,1:N,1:N,ones(1,N));



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
x0 = smat_to_vec(x0_mat)';

t = 1;
MAXITERS = 1000;
niters = 1;
primal_obj = 0;

% Generate posynomial expressions for f0

E = ones(1,size(A,1));

% INITIALIZE AND CHECK FEASIBILITY OF x0

x = x0;  y = A*x;  f = log(E*exp(y))-dot(c,x);


ntiters = 0;
MAXITERS = 250;  % maximum total number of Newton iterations
NTTOL = 1e-3;    % stop centering when Nt decrement <= NTTOL 
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

% Generate equality constraint which says matrix elements should
% sum to zero
A_lineq = lineq_vec;
n_lineq = size(A_lineq,1);

% Correct this
m = NVAR;

while (big_iters <= MAXBIGITERS)

  loglik = 0;
  fprintf(1,'t=%g  N=%4d f=%g. Lik=%g\n',t,big_iters,full(f),full(train_lik));   
   
   center_iters = 0;
   while (center_iters <= MAXITERS)

      % FORM GRADIENT AND HESSIAN

      [bar_grad,bar_hess] = det_barrier_grad(x,F);

      % Generate function for the exponent
%      gradf = exp(y)./(E'*(E*exp(y)));
      % Generate a vector that looks like exp(y)/sum(exp(y))
      mxy = max(y);
      % Start out with the log
      norm_expy = y-mxy-log(sum(exp(y-mxy)));
      norm_expy = exp(norm_expy);
      gradf = norm_expy;
%      gradf = exp(y)./sum(exp(y));
      
      xmat = vec_to_smat(x,F);     

      gradphi = A' * (gradf*t)-t*d_grad +bar_grad + t*trace_grad; 
%      hessphi =  - (A' * (exp(y) * t./sum(exp(y)).^2)) * (exp(y') * A)  ... 
%	  +  (A'*sparse(diag(gradf.*(E'* t)))) * A + bar_hess; %
      hessphi =  - (A' * (norm_expy * t)) * (norm_expy' * A)  ... 
	  +  (A'*sparse(diag(gradf.*(E'* t)))) * A + bar_hess; %
      
      if (rcond(full(hessphi))<eps)
	fprintf(1,'Ill conditioned in Newton\n');
	break;
      end
	
      % Compute the newton direction (non equality constrained)
      

      
      % Compute it for the linear equality constrained case
      dx = -hessphi\gradphi;

      center_iters = center_iters + 1;
      fprime = gradphi'*dx;
      dy = A*dx;
      if (-fprime/2 < NTTOL), break; end; 


      % BACKTRACKING LINE SEARCH


%      phi = t*f - log(det(xmat)); %-log(sum(xmat));
      %%sum(log(-f(2:m+1)));
      mxy = max(y);
      lse = log(sum(exp(y-mxy)))+mxy;
      oldf = lse-dot(c,x) + tr_wgt*trace(xmat);						    
      oldphi = t*oldf - log(det(xmat));
      
      s = 1e-1;
      for lsiters = 1:MAXLSITERS
	newmat = full(vec_to_smat(x+s*dx,F));
	
	% Calculate log(sum(exp(y+s*dy)));
	newy = y+s*dy;
	mx_newy = max(newy);
	lse = log(sum(exp(newy-mx_newy)))+mx_newy;
	newf = lse-dot(c,(x+s*dx)) + tr_wgt*trace(newmat);

         if (min(eig(newmat))> 0) % & sum(newmat)>0
            newphi = t*newf - log(det(newmat)); 
            if (newphi < oldphi + ALPHA*s*fprime)
%            if (newphi < phi )
              break;
            else
%	      newphi-oldphi
	    end
         end
         s = BETA*s;
      end

      x = x+s*dx;
      y = A*x;
      newmat = vec_to_smat(x+s*dx,F);
      mxy = max(y);
      lse = log(sum(exp(y-mxy)))+mxy;
      train_lik = lse-dot(c,x);
      f = train_lik+tr_wgt*trace(newmat);
%      fprintf('The trace is %g\nEigs:',trace(newmat));
%      eig(newmat)
      % Calculate Newton decrement

%      if (abs(f-last_f)<CENTER_MINDELTA)
%          fprintf(1,'Break Newton at %d\n',center_iters);
%          break;
%      end
      last_f = f;
      
     fprintf('In Newton %g. Iter=%d\n',newphi,lsiters);
   end

   t = min(MU*t, (m+1)/ABSTOL);  
   if t>100000
       break;
   end
   big_iters = big_iters+1;
end

% This is (n-1) SDP
m = vec_to_smat(x,F);
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

function v = trans_ind(N,xi,yi)

v = sum(N:-1:N-xi+2)+(yi-xi+1);


