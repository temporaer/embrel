function pxy = equimarg(pxy0)

pxy = pxy0;

while 1
    pxy = make_cond_dist(pxy,0);
    pxy = pxy/sum(pxy(:));
    pxy = make_cond_dist(pxy,1);
    pxy = pxy/sum(pxy(:));
    px = sum(pxy,2);
    st = std(px);
    fprintf(1,'Std(px)=%g\n',st);
    if st<1e-8
        return;
    end
end
    