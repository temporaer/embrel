function [f,g] = code_sdp_grad(x,params)

global const; % = params.const;
global A; % = params.A;
global px0; % = params.px0;
global py0; % = params.py0;
global xinds; % = params.xinds;
global yinds; % = params.yinds;
global NX;
global NY;
t = params.t;
global c; % = params.c;
global trace_grad;
global tr_wgt;

y = A*x+const;

gmat = vec_to_smat_fast(x(1:end-1),NX+NY-1);

bar_grad = det_barrier_grad_gd(gmat);

bar_grad(end+1,1)=0;


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
big_lse = lse;
big_gradf = gradf;



bigvecx = sparse(size(A,1),1);

% Generate gradient and Hessians for px0 constraints
for xi=1:NX
  yx = y(xinds(xi,:))-log(px0(xi));
  mxy = max(yx);
  % Start out with the log
  lse = log(sum(exp(yx-mxy)))+mxy;      
  norm_expy = yx-lse;
  norm_expy = exp(norm_expy);
  gradf = norm_expy;
  px_lse(xi) = lse;
  bigvecx(xinds(xi,:)) = 1/lse*gradf;
end




bigvecy = sparse(size(A,1),1);

for yi=1:NY
  yy = y(yinds(yi,:))-log(py0(yi));
  mxy = max(yy);
  lse = log(sum(exp(yy-mxy)))+mxy;      
  % Start out with the log
  norm_expy = yy-lse; 
  norm_expy = exp(norm_expy);
  gradf = norm_expy;
  py_lse(yi)= lse;
  bigvecy(yinds(yi,:)) = 1/lse*gradf; 
end



[ev,psdness] = chol(gmat);
ev = diag(ev);

%eg = eig(gmat);

%if isinf(abs(min(eg))) | isinf(abs(big_lse)) | ...
%      isinf(abs(max(px_lse))) | isinf(abs(max(py_lse)))
%  disp('Errorororr');
%  keyboard
%end

%if (min(eg)<= 0 | big_lse>=0 | max(px_lse)>=0 | ...
%      max(py_lse)>=0)
if (psdness>0) | big_lse>=0 | max(px_lse)>=0 | max(py_lse)>=0
  f = Inf;
  g=Inf*ones(size(x));
  return;
end

logdet = 2*sum(log(ev));

oldf = dot(c,x) + tr_wgt*trace(gmat);						    
f = t*oldf - logdet -log(-lse)-sum(log(-px_lse))-sum(log(-py_lse));


if nargout>1
  gradphi = t*c + t*trace_grad - 1/big_lse*A' * big_gradf + bar_grad; 
  gradphi = gradphi - A'*bigvecx;
  gradphi = gradphi - A'*bigvecy;
  g = gradphi;
end


