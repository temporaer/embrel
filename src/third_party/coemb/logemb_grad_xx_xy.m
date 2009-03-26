function [f,g] = logemb_grad(x,params)

dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;
pxx0 = params.pxx0;
lambda = params.lambda;

[phi,psi] = read_legrad_x(x,params);
phi_tr = phi';
psi_tr = psi';

% Generate model distribution
dxy = my_dist(phi,psi_tr);

pxy = exp(-dxy);
Zxy = sum(pxy(:));
pxy = pxy/Zxy;

dxx = my_dist(phi,phi_tr);
pxx = exp(-dxx);
Zxx = sum(pxx(:));
pxx = pxx/Zxx;

px = sum(pxy,2);
py = sum(pxy,1);
px0 = sum(pxy0,2);
py0 = sum(pxy0,1);
px_xx = sum(pxx,2);
px_x0 = sum(pxx0,2);

gphi = diag(px)*phi-pxy*psi - (diag(px0)*phi-pxy0*psi);
gphi = gphi+lambda*(diag(px_xx)*phi-pxx*phi - (diag(px_x0)*phi-pxx0*phi));

%gpsi2 = diag(py)*psi-pxy'*phi - (diag(py0)*psi-pxy0'*phi);
gpsi = psi_tr*diag(py)-phi_tr*pxy - (psi_tr*diag(py0)-phi_tr*pxy0);
gpsi = gpsi';

g = [-gphi(:);-gpsi(:)];

f = log(Zxy)+sum(pxy0(:).*dxy(:));
f = f + lambda*(log(Zxx)+sum(pxx0(:).*dxx(:)));















if 0

psi_exp_d = pxy0*psi;
phi_exp_d = phi_tr*pxy0;

for xi=1:NX
    for di=1:dim
        precon_phi(xi,di) = px(xi)-px0(xi);
        tmp = 0;
        for yi=1:NY
            dlogpxy = px(xi)*phi(xi,di)-psi_exp_d(xi,di)-(phi(xi,di)-psi(yi,di));
            tmp = tmp+ pxy(xi,yi)*dlogpxy*(phi(xi,di)-psi(yi,di));
        end
        precon_phi(xi,di) = precon_phi(xi,di)+tmp;
    end
end

for yi=1:NY
    for di=1:dim
        precon_psi(yi,di) = py(yi)-py0(yi);
        tmp = 0;
        for xi=1:NX
            dlogpxy = py(yi)*psi(yi,di)-phi_exp_d(di,yi)-(psi(yi,di)-phi(xi,di));
            tmp = tmp+ pxy(xi,yi)*dlogpxy*(psi(yi,di)-phi(xi,di));
        end
        precon_psi(yi,di) = precon_psi(yi,di)+tmp;
    end
end

gphi = gphi./precon_phi;
gpsi = gpsi./precon_psi;
end
