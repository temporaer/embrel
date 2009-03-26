function  cm = sp_make_cond_dist(m)

NX  = size(m,1);


s = spdiags(1./sum(m,2),0,NX,NX);

cm = s*m;