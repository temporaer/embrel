function [m,train_lik] = code_sdp_solve5(nxy,NX,NY,tr_wgt,bfgs_params)  %le_solve(samp,NX,NY)

force_marg = 1;
% tr_wgt is the weight of the trace constraint

nxy = nxy/sum(nxy(:));

px0 = sum(nxy,2);
py0 = sum(nxy,1);

%F={};
% The F matrices take care of the distance matrix being PSD
N = NX+NY;
NM1 = N-1;
NVAR = (N-1)*N/2;
NDISTS = NX*NY;
trace_grad =sparse(NVAR,1);

% Initialize x0 to NX+NY points radnomly spaced around the origin
rand('seed',0);
[xi,yi,v] = find(nxy);
% Get gradient of the \sum n(x,y)d(x,y) element
%[c] = get_kappa_op_inds(N,xi,NX+yi,v);
%c = c';
% This element is for +log(Z)
%c(end+1) = 1;
% Find how all the dxy are generated from p
% Find how all the dxy are generated from p
[i,j] = ndgrid(1:NX,1:NY);
i = i(:); j=j(:);
A = get_kappa_op_inds(N,i,NX+j,ones(1,length(i)));

[nxyinds,tmp,v] = find(nxy(:));
% Translate to ot
c = v'*A(nxyinds,:);
c = c';
c(end+1:end+NX) = px0;

%[i,j] = meshgrid(1:NX,1:NY);
%i = i(:); j=j(:);
%[tmp,A] = get_kappa_op_inds(N,i,NX+j,ones(1,length(i)));
for xi=1:NX
  xinds(xi,:) = find(i==xi)';
end
for yi=1:NY
  yinds(yi,:) = find(j==yi)';
end

if force_marg
  for k=1:length(i)
    const(k,1) = log(py0(j(k)));
  end
else
  const = 0;
end

A = -A;
% Get gradient of the trace element
%trace_grad = get_kappa_op_inds(N,1:N,1:N,ones(1,N));



%F = sparse((N-1)^2,NVAR);
%F = cell(NVAR,1);
trind = 0;

%lineq_vec = zeros(1,NVAR);

for i=1:N-1
  for j=i:N-1
    trind = trind+1;
    mat = sparse(N-1,N-1);
    mat(i,j) = 1;
    mat(j,i) = 1;
    %        F(:,trind) = mat(:);
%    F{trind} = mat;
    if (j==i)
      trace_grad(trind)=tr_wgt;
    end
  end
end

trace_grad(end+1:end+NX,1)=0;

points = rand(N,NX+NY);
mn = mean(points,1);
points = points - repmat(mn,NX+NY,1);
% Generate their Distance matrix
x0_mat = my_dist(points,points');
clear points;
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
MAXITERS = 1000;
niters = 1;
primal_obj = 0;

% Generate posynomial expressions for f0

%E = ones(1,size(A,1));

% INITIALIZE AND CHECK FEASIBILITY OF x0

% This is for the condition \sum_{x,y} exp(-d(xy)-a) <=1
A(:,end+1:end+NX) = -1*repmat(eye(NX),NY,1);

y = A*x+const;  f =dot(c,x);


ntiters = 0;
MAXITERS = 500;  % maximum total number of Newton iterations
NTTOL = 1e-5;    % stop centering when Nt decrement <= NTTOL 
ALPHA = 0.01; 
BETA = 0.5;      % parameters for backtracking line search
MAXLSITERS = 200; % maxiters in line searches backtracking
MU = 10; %2.5;         % mu for update tplus = mu*t;
ABSTOL=1e-7;     % duality gap based stopping criterion
big_iters = 0;
MAXBIGITERS = 25;
CENTER_MINDELTA = 1e-5;
DIFFTOL = 1e-5;
last_f = 0;
train_lik = 0;

% Generate equality constraint which says matrix elements should
% sum to zero
%A_lineq = lineq_vec;
%n_lineq = size(A_lineq,1);

% Correct this
m = NVAR;
NVAR = NVAR+1;
last_train_lik = 1e9;

while (big_iters <= MAXBIGITERS)

  loglik = 0;
  fprintf(1,'t=%g  N=%4d f=%g. Lik=%g\n',t,big_iters,full(f),full(train_lik));   
   
   center_iters = 0;
   while (center_iters <= MAXITERS)

      % FORM GRADIENT AND HESSIAN

      [bar_grad,bar_hess] = det_barrier_grad(x(1:end-NX),N-1);

      bar_grad(end+1:end+NX,1)=0;
      bar_hess(:,end+1:end+NX) = 0;
      bar_hess(end+1:end+NX,:) = 0;      
      
      % Generate function for the exponent
      % Generate a vector that looks like exp(y)/sum(exp(y))
      mxy = max(y);
      % Start out with the log
      norm_expy = y-mxy-log(sum(exp(y-mxy)));
      norm_expy = exp(norm_expy);
      lse = log(sum(exp(y-mxy)))+mxy;      
      gradf = norm_expy;
%      gradf = exp(y)./sum(exp(y));
      
      xmat = vec_to_smat_fast(x(1:end-NX),N-1);     

%      gradphi = t*c + t*trace_grad - 1/lse*A' * gradf + bar_grad; 
      gradphi = t*c + t*trace_grad + bar_grad ; 
      % Next is the Hessian which is non-zero only for the barrier functions

      mdl = 1/lse*norm_expy * norm_expy'+1/lse^2*gradf*gradf'-1/lse*spdiags(gradf,0,NDISTS,NDISTS);
%      hessphi = A'*mdl*A + bar_hess;
      hessphi =  bar_hess;
      

      
      if force_marg
	
	% Generate gradient and Hessians for px0 constraints
	for xi=1:NX
	  Ax = A(xinds(xi,:),:);
	  yx = y(xinds(xi,:));
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
      
      [ev,posdef] = chol(xmat);
      ev = diag(ev);
      logdet = 2*sum(log(ev));						    
      
      if force_marg
	oldphi = t*oldf - logdet-sum(log(-px_lse));
	old_pxlse = px_lse;
      else
	oldphi = t*oldf - logdet-log(-lse);
      end
      

      s = 1e-1;
      for lsiters = 1:MAXLSITERS
	newmat = full(vec_to_smat_fast(x(1:end-NX+1)+s*dx(1:end-NX+1),N-1));
	
	% Calculate log(sum(exp(y+s*dy)));
	newy = y+s*dy;
	newx = x+s*dx;
	dist_m = kappav(newmat);

        newlogz = newx(end-NX+1:end);
        pxy = exp(-dist_m(1:NX,NX+1:end)-repmat(newlogz,1,NY)+repmat(log(py0),NX,1));	
	
	mx_newy = max(newy);
	lse = log(sum(exp(newy-mx_newy)))+mx_newy;
	if force_marg
	  for xi=1:NX
	    curr_y = newy(xinds(xi,:));
	    mx_curry = max(curr_y);	  
	    px_lse(xi) = log(sum(exp(curr_y-mx_curry)))+mx_curry;	  
	  end
	end
	newf = dot(c,(x+s*dx)) + tr_wgt*trace(newmat);
	
	[ev,psdness] = chol(newmat);
	ev = diag(ev);
	logdet = 2*sum(log(ev));
	
	if (psdness==0 &  (max(px_lse)<0))
	  newphi = t*newf - logdet  -sum(log(-px_lse));
	   if (newphi < oldphi + ALPHA*s*fprime)
              break;
            end
         end
         s = BETA*s;
      end

      x = x+s*dx;
      y = A*x+const;
      newmat = vec_to_smat_fast(x(1:end-NX+1),N-1);
      mxy = max(y);
      lse = log(sum(exp(y-mxy)))+mxy;
      train_lik = dot(c,x);
      f = train_lik+tr_wgt*trace(newmat);
      fprintf('In Newton %g. Iter=%d\n',newphi,center_iters);
   end

   t = min(MU*t, (m+1)/ABSTOL);  
   rt =   abs((last_train_lik-train_lik)/last_train_lik);
   fprintf('Rat=%g\n',rt);
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

%for k=0:length(x)-1
%    v = [v diag(x,k)];
%end

logzx = x(end-NX+1:end);
pxy = exp(-dist_m(1:NX,NX+1:end)-repmat(logzx,1,NY)+repmat(log(py0),NX,1));


y = pxy;
% Start out with the log
lse = full(sum(y,2));      

fprintf('NormTo=%g-%g\n',min(lse),max(lse));


function v = trans_ind(N,xi,yi)

v = sum(N:-1:N-xi+2)+(yi-xi+1);


