function [f,g] = precond_code_sdp_grad(x,params)

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

x0 = x;
gmat = vec_to_smat_fast(x(1:end-1),NX+NY-1);

dmat = kappav(gmat);
dmat = dmat(1:NX,NX+1:end);
dmat = -dmat +log(px0)*ones(1,NY)+ones(NX,1)*log(py0)-x(end);

y = dmat(:);




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



bigvecx = zeros(length(big_gradf),1);

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




bigvecy = zeros(length(big_gradf),1);

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


if  isinf(abs(big_lse)) | ...
      isinf(abs(max(px_lse))) | isinf(abs(max(py_lse)))
  disp('Errorororr');
  keyboard
end

if (psdness>0) | big_lse>=0 | max(px_lse)>=0 | max(py_lse)>=0
  f = Inf;
  g=Inf*ones(size(x));
  return;
end

logdet = 2*sum(log(ev));

oldf = dot(c,x) + tr_wgt*trace(gmat);						    
f = t*oldf - logdet -log(-lse)-sum(log(-px_lse))-sum(log(-py_lse));


if nargout==1
   return;
end



%other = A'*big_gradf;
%other = vec_to_smat_fast(other,N);
exA = zeros(N,N);
d3 = reshape(big_normexpy,NX,NY);

n=N+1;
x = -1/(n+sqrt(n));
y = -1/sqrt(n);
a_diff = 2*(x-y)^2;
b_diff = 2*(x+1-y)*(x-y);
a_eq = (x-y)^2;
b_eq = (x+1-y)^2;




exA(1:NX-1,NX:end) = -2*d3(2:end,:);
add = ones(NX-1,NY)*sum(d3(1,:))*a_diff+repmat(d3(1,:),NX-1,1)*2*(x-y);
exA(1:NX-1,NX:end) = exA(1:NX-1,NX:end) + add;

exA(1:NX-1,1:NX-1) = sum(d3(1,:))*a_diff;

constx = (a_eq-a_diff)*sum(d3(1,:));
consty = a_eq*sum(d3(1,:));

diagX = sum(d3(2:end,:),2)+constx;

exA = exA + diag([diagX' zeros(1,NY)]);

exA(NX:end,NX:end) = sum(d3(1,:))*a_diff;
diagY = sum(d3(2:end,:),1)+constx;
diagY = diagY+(b_eq - a_eq)*d3(1,:);
exA = exA + diag([zeros(1,NX-1) diagY]);

for yi1=NX:N
    for yi2 = yi1+1:N
       exA(yi1,yi2) = exA(yi1,yi2) + (d3(1,yi1-NX+1)+d3(1,yi2-NX+1))*2*(x-y);
   end
end

exB = zeros(N,N);
dx = reshape(bigvecx,NX,NY);
%otherX = vec_to_smat_fast(A'*bigvecx,N);

exB(1:NX-1,NX:end) = -2*dx(2:end,:);
add = ones(NX-1,NY)*sum(dx(1,:))*a_diff+repmat(dx(1,:),NX-1,1)*2*(x-y);
exB(1:NX-1,NX:end) = exB(1:NX-1,NX:end) + add;

exB(1:NX-1,1:NX-1) = sum(dx(1,:))*a_diff;

constx = (a_eq-a_diff)*sum(dx(1,:));
consty = a_eq*sum(dx(1,:));

diagX = sum(dx(2:end,:),2)+constx;

exB = exB + diag([diagX' zeros(1,NY)]);

exB(NX:end,NX:end) = sum(dx(1,:))*a_diff;
diagY = sum(dx(2:end,:),1)+constx;
diagY = diagY+(b_eq - a_eq)*dx(1,:);
exB = exB + diag([zeros(1,NX-1) diagY]);

for yi1=NX:N
    for yi2 = yi1+1:N
       exB(yi1,yi2) = exB(yi1,yi2) + (dx(1,yi1-NX+1)+dx(1,yi2-NX+1))*2*(x-y);
   end
end


exC = zeros(N,N);
dy = reshape(bigvecy,NX,NY);
%otherY = vec_to_smat_fast(A'*bigvecy,N);

exC(1:NX-1,NX:end) = -2*dy(2:end,:);
add = ones(NX-1,NY)*sum(dy(1,:))*a_diff+repmat(dy(1,:),NX-1,1)*2*(x-y);
exC(1:NX-1,NX:end) = exC(1:NX-1,NX:end) + add;

exC(1:NX-1,1:NX-1) = sum(dy(1,:))*a_diff;

constx = (a_eq-a_diff)*sum(dy(1,:));
consty = a_eq*sum(dy(1,:));

diagX = sum(dy(2:end,:),2)+constx;

exC = exC + diag([diagX' zeros(1,NY)]);

exC(NX:end,NX:end) = sum(dy(1,:))*a_diff;
diagY = sum(dy(2:end,:),1)+constx;
diagY = diagY+(b_eq - a_eq)*dy(1,:);
exC = exC + diag([zeros(1,NX-1) diagY]);

for yi1=NX:N
    for yi2 = yi1+1:N
       exC(yi1,yi2) = exC(yi1,yi2) + (dy(1,yi1-NX+1)+dy(1,yi2-NX+1))*2*(x-y);
   end
end


bar_grad = det_barrier_grad_gd(gmat);
bar_grad(end+1,1)=0;

exA = [smat_to_vec_c(exA') 1]';
exB = [smat_to_vec_c(exB') sum(bigvecx)]';
exC = [smat_to_vec_c(exC') sum(bigvecy)]';

gradphi = gmat*vec_to_smat_fast(t*c(1:end-1) + t*trace_grad(1:end-1) + 1/big_lse*exA(1:end-1)+exB(1:end-1)+exC(1:end-1),length(gmat)) + t*eye(length(gmat)); 

g = smat_to_vec_c(gradphi)';
g(end+1) = 1/big_lse*exA(end)+exB(end)+exC(end)+t*c(end)+t*trace_grad(end);

gradphi2 = t*c + t*trace_grad + 1/big_lse*exA - bar_grad; 
gradphi2 = gradphi2 + exB;
gradphi2 = gradphi2 + exC;

v = vec_to_smat_fast(gradphi2(1:end-1),length(gmat));
% Multiply off diagonals by half, to make this the derivative w.r.t
% the full Gram matrix
v = gmat*v;
keyboard
v = smat_to_vec_c(v);
g = [v gradphi2(end)]';

%gradphi2(1:end-1,1:end-1) = x*gradphi2(1:end-1,1:end-1);
%keyboard

