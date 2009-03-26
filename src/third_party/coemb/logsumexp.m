function lse = logsumexp(logf)

base = 1+1e-5;
mxlogf = max(logf);  

lse = log(sum(sum(base.^(logpxy-mxy))));
