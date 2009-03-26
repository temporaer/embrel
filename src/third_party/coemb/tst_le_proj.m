%function tst_le_proj(nxy)

% Compare projection solution to direct

dims = [2:2:20];
reps = 3;

for di=1:length(dims)
    for ri=1:reps
        [phi,psi] = run_le(nxy,dims(di));
        lk_direct(ri,di) = get_le_lik(phi,psi,nxy);
%        [phi,psi] = run_le_proj(nxy,dims(di));
%        lk_proj(ri,di) = get_le_lik(phi,psi,nxy);
        [phi,psi] = run_sdp_proj(nxy,dims(di));
        lk_sdp(ri,di) = get_le_lik(phi,psi,nxy);
    end
end

errorbar(dims,mean(lk_direct),std(lk_direct,[],1));
%hold on
%ccerrorbar(dims,mean(lk_proj),std(lk_proj,[],1),'g');
%hold on
%errorbar(dims,mean(lk_sdp),std(lk_sdp,[],1),'r');
