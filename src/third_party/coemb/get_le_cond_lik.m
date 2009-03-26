function  [lik,pxy] = get_le_cond_lik(phi,psi,nxy,b_ycond)

nxy = nxy/sum(nxy(:));
py0 = sum(nxy,1);
px0 = sum(nxy,2);

NX = size(nxy,1);
NY = size(nxy,2);

t = [phi;psi]*[phi;psi]';
d = g_to_d(t);
d = my_dist(phi,psi');
subd = d; %(1:NX,NX+1:end);
%found_pxy = exp(-subd);
%found_pxy = found_pxy(1:NX,NX+1:end);
%found_pxy = found_pxy/sum(found_pxy(:));
%logpxy = log(found_pxy); %-curr_dmat-log(Z);  %pxy/sum(pxy(:));
px0 = sum(nxy,2);
py0 = sum(nxy,1);

mx = mean(subd(:));
if b_ycond
  logpxy = repmat(log(px0),1,NY) + mx-subd-log(sum(sum(exp(-subd+mx+repmat(log(px0),1,NY)))));  
  pxy = exp(logpxy);
  pxy = make_cond_dist(pxy,1);
else
  logpxy =  repmat(log(py0),NX,1) + mx-subd-log(sum(sum(exp(-subd+mx+repmat(log(py0),NX,1)))));  
  pxy = exp(logpxy);
  pxy = pxy+rand(size(pxy))*eps;
  pxy = make_cond_dist(pxy,0);
end

logpxy = log(pxy);

lik = sum(sum(nxy.*logpxy));

