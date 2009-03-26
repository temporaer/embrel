function [x] = logemb_proj(x,params)

dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;

%[v,j] = max(pxy0(:));
%[i,j] = ind2sub(size(pxy0),j);

[phi,psi] = read_legrad_x(x,params);
%phi = phi-repmat(mean(phi),NX,1);
%psi = psi-repmat(mean(psi),NY,1);

%a=1;
%phi = phi*diag(sqrt(a)./sqrt(sum(phi.^2)));
%psi = psi*diag(sqrt(a)./sqrt(sum(psi.^2)));

if 1
mn = mean([phi;psi]);
phi=phi-repmat(mn,NX,1);
psi=psi-repmat(mn,NY,1);

%sq = sqrt(sum([phi;psi].^2));
%min(sq)
%sq = sq*0+1; %min(sq);
%phi = phi*diag(1./sq);
%psi = psi*diag(1./sq);
%v = sqrt(sum([phi(:);psi(:)].^2));
v = std([phi(:);psi(:)])/0.1;
phi = phi/v;
psi = psi/v;
end

x = [phi(:);psi(:)];