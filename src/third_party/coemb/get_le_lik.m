function  [lik,pxy] = get_le_lik(phi,psi,nxy,do_marg)

if nargin<4
  do_marg = 1;
end

nxy = nxy/sum(nxy(:));

NX = size(nxy,1);

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
%logZ = log(sum(exp(-subd(:)+mx)))-mx;
if do_marg
  logpxy = mx-subd-log(sum(sum((px0*py0).*exp(-subd+mx))))+log(px0*py0);
else
  logpxy = mx-subd-log(sum(sum(exp(-subd+mx))));  
end
lik = sum(sum(nxy.*logpxy));

pxy = exp(logpxy);