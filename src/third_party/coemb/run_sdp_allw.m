function [phi,psi,liks,evs] = run_sdp_allw(nxy,curr_dim,tr_wgts)

% Run in full dimension and project
nxy = nxy/sum(nxy(:));

for wi=1:length(tr_wgts)
    w = tr_wgts(wi);
    [phis{wi},psis{wi},ev] = run_sdp_proj(nxy,curr_dim,tr_wgts(wi));
    evs(:,wi) = ev;
    liks(wi) = get_le_lik(phis{wi},psis{wi},nxy);
    fprintf(1,'Weight=%g Lik=%g\n',w,liks(wi));
end    

% Find best solution
[v,i] = max(liks);
phi = phis{i};
psi = psis{i};