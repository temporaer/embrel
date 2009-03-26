function [f,g] = nomarg_sdp_grad(x,params)

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

gmat = vec_to_smat_fast(x(1:end-1),NX+NY-1);

%disp('B');

dmat = kappav(gmat);
dmat = dmat(1:NX,NX+1:end);
dmat = -dmat -x(end);

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

%disp('D');

[ev,psdness] = chol(gmat);
ev = diag(ev);

%disp('G');

if  isinf(abs(big_lse))
  disp('Errorororr');
end

if (psdness>0) | big_lse>=0
  f = Inf;
  g=Inf*ones(size(x));
  return;
end

%disp('H');

logdet = 2*sum(log(ev));

oldf = dot(c,x) + tr_wgt*trace(gmat);						    
f = t*oldf - logdet -log(-lse);


if nargout==1
   return;
end


%disp('I');

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

%disp('J');


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

exA = [smat_to_vec_c(exA') 1]';
exA = exA/big_lse;

%disp('L4');

bar_grad = det_barrier_grad_gd(gmat);
bar_grad(end+1,1)=0;

%disp('L5');

gradphi = t*c + t*trace_grad + exA + bar_grad; 
g = gradphi;

%disp('M');
