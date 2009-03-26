function [phi,psi,lik] = sdp_proj_optweights(nxy)

wts = [0 0.001:0.0005:0.01]; %1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2
                                                                             % 1e-1 1]; % 2 5 10 100 1000];
%wts = [1e-4:1e-4:1e-3];
wts = [0 0.001 0.01 0.1 1 10 100];
wts = logspace(-4,2,20);
%wts = [0 0.01]; 									     
%[phi,psi] = run_le_proj(nxy,2);
%base_lik = get_le_lik(phi,psi,nxy);

%nxy = pxy0;
%clear pxy0;

for wi = 1:length(wts)
    w = wts(wi);
    [phi_w{wi},psi_w{wi}] = run_sdp_proj(nxy,2,w);
    rep= get_le_lik(phi,psi,nxy);
    sdp_liks(wi) = max(rep);
end

[lik,bst_ind] = max(rep);
phi = phi_w{wi};
psi = psi_w{wi};
