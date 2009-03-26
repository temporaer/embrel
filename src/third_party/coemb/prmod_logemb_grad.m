function [f,g] = prmod_logemb_grad(x,params)

dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;
%k = params.k;
bet = params.bet;
sig = params.sig;

k = 1/(2+1/(bet*sig^2));
%[v,j] = max(pxy0(:));
%[i,j] = ind2sub(size(pxy0),j);

[phi,psi] = read_legrad_x(x,params);
phi_tr = phi';
psi_tr = psi';

% Generate model distribution
d = my_dist(phi,psi_tr);

dd = bet*((1-k)/2*d+(1-2*k)*phi*psi_tr);
mx = max(dd(:));
logpxy = mx-dd-log(sum(sum(exp(-dd(:)+mx))));
pxy = exp(logpxy);

%pxy = exp(-bet*((1-k)/2*d-(1-2*k)*phi*psi_tr));

%Z = sum(pxy(:));
%pxy = pxy/Z;

px = sum(pxy,2);
py = sum(pxy,1);
px0 = sum(pxy0,2);
py0 = sum(pxy0,1);

gphi = (1-k)*diag(px)*phi-k*pxy*psi - ((1-k)*diag(px0)*phi-k*pxy0*psi);
%gpsi2 = diag(py)*psi-pxy'*phi - (diag(py0)*psi-pxy0'*phi);
gpsi = (1-k)*psi_tr*diag(py)-k*phi_tr*pxy - ((1-k)*psi_tr*diag(py0)-k*phi_tr*pxy0);
gpsi = gpsi';

g = [-gphi(:);-gpsi(:)];

f = -sum(pxy0(:).*log(pxy(:)));
%log(Z)+sum(pxy0(:).*d(:));















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
