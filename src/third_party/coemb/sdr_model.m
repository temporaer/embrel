function p = sdr_model(phi,psi,a,b,bcond)

NY = length(b);
NX = length(a);

p = exp(-my_dist(phi,psi')+a*ones(1,NY)+ones(NX,1)*b);
if bcond
  p = make_cond_dist(p,0);
else
  p = p/sum(p(:));
end