function [f,g] = cond_sdp_grad(x,params)

global const; % = params.const;
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
global pxy0;

%y = A*x+const;

N = NX+NY-1;

%disp('A');

gmat = vec_to_smat_fast(x(1:end-NX),NX+NY-1);

%disp('B');

dmat = kappav(gmat);
dmat = dmat(1:NX,NX+1:end);
logzx = x(end-NX+1:end);
dmat = -dmat - repmat(logzx,1,NY)+repmat(log(py0),NX,1);

y = dmat(:);


%disp('C');

% Generate function for the exponent
%      gradf = exp(y)./(E'*(E*exp(y)));
% Generate a vector that looks like exp(y)/sum(exp(y))
mxy = max(y);
% Start out with the log
norm_expy = y-mxy-log(sum(exp(y-mxy)));
norm_expy = exp(norm_expy);
big_normexpy = norm_expy;
lse = log(sum(exp(y-mxy)))+mxy;      
gradf = norm_expy;
%      gradf = exp(y)./sum(exp(y));
big_lse = lse;
big_gradf = gradf;

clear gradf 
%disp('D');

bigvecx = zeros(length(big_gradf),1);

% Generate gradient and Hessians for px0 constraints
for xi=1:NX
  yx = y(xinds(xi,:));
  mxy = max(yx);
  % Start out with the log
  lse = log(sum(exp(yx-mxy)))+mxy;      
  norm_expy = yx-lse;
  norm_expy = exp(norm_expy);
  gradf = norm_expy;
  px_lse(xi) = lse;
  bigvecx(xinds(xi,:)) = 1/lse*gradf;
end

clear gradf norm_expy yx 

%disp('E');

%disp('F');

[ev,psdness] = chol(gmat);
ev = diag(ev);

%disp('G');

if  isinf(abs(big_lse)) | isinf(abs(max(px_lse)))
  disp('Errorororr');
end

if (psdness>0) | max(px_lse)>=0 
  f = Inf;
  g=Inf*ones(size(x));
  return;
end

%disp('H');

logdet = 2*sum(log(ev));

oldf = dot(c,x) + tr_wgt*trace(gmat);						    
f = t*oldf - logdet -sum(log(-px_lse));


if nargout==1
   return;
end


%disp('I');

dx = reshape(bigvecx,NX,NY);
exA = make_c(dx,NX,NY,N);
exA = [smat_to_vec_c(exA') 1./px_lse]';


bar_grad = det_barrier_grad_gd(gmat);
bar_grad(end+1:end+NX,1)=0;

%disp('L5');

gradphi = t*c + t*trace_grad + exA + bar_grad; 
g = gradphi;

%disp('M');
