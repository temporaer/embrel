function [phi,psi] = norm_features(phi,psi)

NX = size(phi,1);
NY = size(psi,1);

%mn = mean([phi;psi]);
%phi=phi-repmat(mn,NX,1);
%psi=psi-repmat(mn,NY,1);

mn = mean(phi);
phi=phi-repmat(mn,NX,1);

mn = mean(psi);
psi=psi-repmat(mn,NY,1);

sq = std(phi);
phi = phi*diag(1./sq);

sq = std(psi);
psi = psi*diag(1./sq);

%sq = sqrt(sum([phi;psi].^2));
%min(sq)
%sq = sq*0+1; %min(sq);
%phi = phi*diag(1./sq);
%psi = psi*diag(1./sq);
%v = sqrt(sum([phi(:);psi(:)].^2));
%v = std([phi(:);psi(:)]);
%phi = phi/v;
%psi = psi/v;
