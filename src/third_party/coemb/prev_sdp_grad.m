function [f,g,basef] = sdp_grad(x,params)

A = params.A;
c = params.c;
trace_grad = params.trace_grad;
F = params.F;
tr_wgt = params.tr_wgt;
t = params.t;

y = A*x;
gradf = exp(y-y(1))./sum(exp(y-y(1)));

bar_grad = det_barrier_grad_gd(x,F);
newmat = F*x;  %vec_to_smat(x,F);
N = sqrt(size(F,1));
newmat = reshape(newmat,N,N);

if (min(eig(newmat))> 0)
    basef = log(sum(exp(y-y(1))))+y(1)-dot(c,x); 
    newf =  basef + tr_wgt*trace(newmat);    
    f = t*newf - log(det(newmat));     
else
    f = Inf
end
    
g = A' * (gradf*t)-t*c +bar_grad + trace_grad; 
