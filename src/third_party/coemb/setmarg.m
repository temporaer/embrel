function pxy = setmarg(pxy0,px0,py0)

pxy = pxy0;

while 1
    pxy = make_cond_dist(pxy,0);
    pxy = diag(px0)*pxy; %*diag(py0);    
%    pxy = pxy/sum(pxy(:));
    pxy = make_cond_dist(pxy,1);
    pxy = pxy*diag(py0);
%    pxy = pxy/sum(pxy(:));
    px = sum(pxy,2);
    py = sum(pxy,1);
%    fprintf(1,'Diff(px)=%g Diff(py)=%g\n',mean(abs(px-px0)),mean(abs(py-py0)));
    st = max(mean(abs(px-px0)),mean(abs(py-py0)));
    if st<1e-8
        return;
    end
end
    