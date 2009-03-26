

wts = [0 1 2 5 10 100 1000];
[phi,psi] = run_le_proj(nxy,2);
base_lik = get_le_lik(phi,psi,nxy);

for wi = 1:length(wts)
    w = wts(wi);
    [phi,psi] = run_sdp_proj(nxy,2,w);
    sdp_liks(wi) = get_le_lik(phi,psi,nxy);
end

plot([1 length(wts)],[base_lik base_lik],1:length(wts),sdp_liks);