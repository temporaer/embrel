function [phi,psi] = norm_sdr_features(phi,psi,px,py)

NX = size(phi,1);
NY = size(psi,1);

%mn = px*phi;
%phi  =phi - repmat(mn,NX,1);

%mn = py*psi;
%psi = psi -repmat(mn,NY,1);

X = phi*psi'.*(sqrt(px)*sqrt(py));
d = size(phi,2);

[U,S,V] = svds(X,d);


phi = diag(1./sqrt(px))*U;
psi = diag(1./sqrt(py))*V;

return;






%mn = mean([phi;psi]);
%phi=phi-repmat(mn,NX,1);
%psi=psi-repmat(mn,NY,1);


sq = sqrt(std(psi)./std(phi));
phi = phi*diag(sq);
psi = psi*diag(1./sq);
%sq = std(psi);
%psi = psi*diag(1./sq);

%sq = sqrt(sum([phi;psi].^2));
%min(sq)
%sq = sq*0+1; %min(sq);
%phi = phi*diag(1./sq);
%psi = psi*diag(1./sq);
%v = sqrt(sum([phi(:);psi(:)].^2));
%v = std([phi(:);psi(:)]);
%phi = phi/v;
%psi = psi/v;
