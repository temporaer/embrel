function  [lik,Z,pxy] = get_le_lik_copula(phi,psi,nxy)

nxy = nxy/sum(nxy(:));

NX = size(nxy,1);
px = sum(nxy,2);
py = sum(nxy,1);

t = [phi;psi]*[phi;psi]';
d = g_to_d(t);
subd = d(1:NX,NX+1:end);
pxy = (px*py).*exp(-subd);
Z = sum(pxy(:));
pxy = pxy/sum(pxy(:));
%found_pxy = exp(-subd);
%found_pxy = found_pxy(1:NX,NX+1:end);
%found_pxy = found_pxy/sum(found_pxy(:));
%logpxy = log(found_pxy); %-curr_dmat-log(Z);  %pxy/sum(pxy(:));
%mx = max(subd(:));
%margprod = px*py;
%logZ = log(sum(exp(-subd(:)+mx)))-mx+;
%logpxy = mx-subd-log(sum(sum(exp(-subd+mx))))+log(px*py);
%lik = sum(sum(nxy.*log(exp(-subd)/Z)));
lik = sum(nxy(:).*log(pxy(:)));