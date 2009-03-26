function [xc] = proj_grad(fgrad,fproj,par,x0,niter)

alpha0 = 1e-3;
jitter0 = 0;
xc = x0;
alphac = alpha0;

min_del = 1e-12;

for it=1:niter
    [fc,gc] = feval(fgrad,xc,par);
    lpar = [16 2 0.1 0.01 5];    
%    [alphac] = linesearch(fgrad,par, xc,fc,gc, -gc, lpar);
    alphac = alpha0/sqrt(it);
    jitterc = jitter0/sqrt(it);
    xc = xc-alphac*gc;
%    xc = xc + jitterc*2*(0.5-rand(size(xc)));
    
    if ~isempty(fproj)
        xc = feval(fproj,xc,par);
    end
    
    fprintf(1,'It=%d Val=%g\n',it,fc);
    val_hist(it) = fc;
    if (it>5 & abs(val_hist(it-1)-val_hist(it))<min_del)
        return;
    end
        
%    if it>1
%        dval = abs(val_hist(it-1)-val_hist(it));
%    end
%    if it>10 dval<min_del*(min_del+norm(xc,inf))
%        return;
%    end
end

1
%clf
%plot(val_hist);
%pause