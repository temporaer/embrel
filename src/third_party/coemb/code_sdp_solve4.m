function [m,train_lik] = code_sdp_solve4(nxy,NX,NY,tr_wgt)  %le_solve(samp,NX,NY)

% tr_wgt is the weight of the trace constraint

nxy = nxy/sum(nxy(:));

px0 = sum(nxy,2);
py0 = sum(nxy,1);

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
[c] = get_kappa_op_inds(N,xi,NX+yi,v);
c = c';
% This element is for +log(Z)
c(end+1) = 1;
% Find how all the dxy are generated from p
[i,j] = meshgrid(1:NX,1:NY);
i = i(:); j=j(:);
[tmp,A] = get_kappa_op_inds(N,i,NX+j,ones(1,length(i)));
for xi=1:NX
  xinds(xi,:) = find(i==xi)';
end
for yi=1:NY
  yinds(yi,:) = find(j==yi)';
end

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
x0 = smat_to_vec(x0_mat)';

% Add a=log(Z)
x0 = [x0; 20];

t = 1;
MAXITERS = 1000;
niters = 1;
primal_obj = 0;

% Generate posynomial expressions for f0

%E = ones(1,size(A,1));

% INITIALIZE AND CHECK FEASIBILITY OF x0

% This is for the condition \sum_{x,y} exp(-d(xy)-a) <=1
A(:,end+1) = -1;

x = x0;  y = A*x;  f =dot(c,x);


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

      [bar_grad,bar_hess] = det_barrier_grad(x(1:end-1),F);

      bar_grad(end+1,1)=0;
      bar_hess(:,end+1) = 0;
      bar_hess(end+1,:) = 0;      
      
      % Generate function for the exponent
%      gradf = exp(y)./(E'*(E*exp(y)));
      % Generate a vector that looks like exp(y)/sum(exp(y))
      mxy = max(y);
      % Start out with the log
      norm_expy = y-mxy-log(sum(exp(y-mxy)));
      norm_expy = exp(norm_expy);
      lse = log(sum(exp(y-mxy)))+mxy;      
      gradf = norm_expy;
%      gradf = exp(y)./sum(exp(y));
      
      xmat = vec_to_smat(x(1:end-1),F);     

      gradphi = t*c + t*trace_grad - 1/lse*A' * gradf + bar_grad; 
      % Next is the Hessian which is non-zero only for the barrier functions
      hessphi =  -1/lse*A'*sparse(diag(gradf))* A + 1/lse* (A' * ...
	  norm_expy * norm_expy' * A)  + 1/lse^2* ...
	  A'*(gradf*gradf')*A + bar_hess;

      % Generate gradient and Hessians for px0 constraints
      for xi=1:NX
	Ax = A(xinds(xi,:),:);
	yx = Ax*x-log(px0(xi));
	mxy = max(yx);
	% Start out with the log
	norm_expy = yx-mxy-log(sum(exp(yx-mxy)));
	norm_expy = exp(norm_expy);
	lse = log(sum(exp(yx-mxy)))+mxy;      
	gradf = norm_expy;
	
	px_lse(xi) = lse;
	gradphi = gradphi - 1/lse*Ax'*gradf;
	hessphi = hessphi -1/lse*Ax'*sparse(diag(gradf))* Ax + 1/lse* (Ax' * ...
	    norm_expy * norm_expy' * Ax)  + 1/lse^2* ...
	    Ax'*(gradf*gradf')*Ax;
      end

      for yi=1:NY
	Ay = A(yinds(yi,:),:);
	yy = Ay*x-log(py0(yi));
	mxy = max(yy);
	% Start out with the log
	norm_expy = yy-mxy-log(sum(exp(yy-mxy)));
	norm_expy = exp(norm_expy);
	lse = log(sum(exp(yy-mxy)))+mxy;      
	gradf = norm_expy;
	
	py_lse(yi)= lse;
	gradphi = gradphi - 1/lse*Ay'*gradf;
	hessphi = hessphi -1/lse*Ay'*sparse(diag(gradf))* Ay + 1/lse* (Ay' * ...
	  norm_expy * norm_expy' * Ay)  + 1/lse^2* ...
	  Ay'*(gradf*gradf')*Ay;
      end
      

      if (rcond(full(hessphi))<eps)
	fprintf(1,'Ill conditioned in Newton\n');
	break;
      end
	
      % Compute the newton direction
      dx = -hessphi\gradphi;

      center_iters = center_iters + 1;
      fprime = gradphi'*dx;
      dy = A*dx;
      if (-fprime/2 < NTTOL), break; end; 


      % BACKTRACKING LINE SEARCH
      mxy = max(y);
      lse = log(sum(exp(y-mxy)))+mxy;
      oldf = dot(c,x) + tr_wgt*trace(xmat);						    
      oldphi = t*oldf - log(det(xmat))-log(-lse)-sum(log(-px_lse))-sum(log(-py_lse));
      old_pxlse = px_lse;
      old_pylse = py_lse;
      
      s = 1e-1;
      for lsiters = 1:MAXLSITERS
	newmat = full(vec_to_smat(x(1:end-1)+s*dx(1:end-1),F));
	
	% Calculate log(sum(exp(y+s*dy)));
	newy = y+s*dy;
	newx = x+s*dx;
	dist_m = kappav(newmat);
	pxy = exp(-dist_m(1:NX,NX+1:end)-newx(end));	
	
	mx_newy = max(newy);
	lse = log(sum(exp(newy-mx_newy)))+mx_newy;
	for xi=1:NX
	  curr_y = newy(xinds(xi,:))-log(px0(xi));
	  mx_curry = max(curr_y);	  
	  px_lse(xi) = log(sum(exp(curr_y-mx_curry)))+mx_curry;	  
	end

	for yi=1:NY
	  curr_y = newy(yinds(yi,:))-log(py0(yi));
	  mx_curry = max(curr_y);	  
	  py_lse(yi) = log(sum(exp(curr_y-mx_curry)))+mx_curry;	  
	end
	
	newf = dot(c,(x+s*dx)) + tr_wgt*trace(newmat);

         if (min(eig(newmat))> 0 & lse<0 & max(px_lse)<0 & max(py_lse)<0) % & sum(newmat)>0
            newphi = t*newf - log(det(newmat)) - log(-lse)-sum(log(-px_lse))-sum(log(-py_lse));; 
            if (newphi < oldphi + ALPHA*s*fprime)
              break;
            end
         end
         s = BETA*s;
      end

      x = x+s*dx;
      y = A*x;
      newmat = vec_to_smat(x(1:end-1)+s*dx(1:end-1),F);
      mxy = max(y);
      lse = log(sum(exp(y-mxy)))+mxy;
      train_lik = dot(c,x);
      f = train_lik+tr_wgt*trace(newmat);
      fprintf('In Newton %g. Iter=%d\n',newphi,lsiters);
   end

   t = min(MU*t, (m+1)/ABSTOL);  
   if t>100000
       break;
   end
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

pxy = exp(-dist_m(1:NX,NX+1:end)-x(end));
y = pxy(:);
mxy = max(y);
% Start out with the log
lse = log(sum(y));      

fprintf('NormTo=%g\n',exp(lse));
fprintf('%s\n',num2str(sum(pxy,1)));
fprintf('%s\n',num2str(py0));
fprintf('%s\n',num2str(sum(pxy,2)'));
fprintf('%s\n',num2str(px0'));


function v = trans_ind(N,xi,yi)

v = sum(N:-1:N-xi+2)+(yi-xi+1);


