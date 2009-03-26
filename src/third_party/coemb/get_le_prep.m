function  lik = get_le_prep(phi,psi,nxy)

nxy = nxy/sum(nxy(:));

NX = size(nxy,1);

t = [phi;psi]*[phi;psi]';
d = g_to_d(t);
subd = d(1:NX,NX+1:end);
%found_pxy = exp(-subd);
%found_pxy = found_pxy(1:NX,NX+1:end);
%found_pxy = found_pxy/sum(found_pxy(:));
%logpxy = log(found_pxy); %-curr_dmat-log(Z);  %pxy/sum(pxy(:));
logpxy = subd(1,1)-subd-log(sum(sum(exp(-subd+subd(1,1)))));
pxy = exp(logpxy);
logpx = log(sum(pxy,2));
logpxy = logpxy - repmat(logpx,1,NX);

lik = sum(sum(nxy.*logpxy));
