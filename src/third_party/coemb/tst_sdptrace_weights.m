

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
if 1
reps = 50;
    
for ri=1:reps
    [phi,psi] = run_le(nxy,2);
    opt_lik(ri) = get_le_lik(phi,psi,nxy);
    fprintf(1,'OptLik = %g\n',full(opt_lik(ri)));
end
end
reps = 1;

for wi = 1:length(wts)
    w = wts(wi);
    for ri=1:reps
        [phi,psi] = run_sdp_proj(nxy,2,w);
         rep(ri)= get_le_lik(phi,psi,nxy);
     end
    sdp_liks(wi) = max(rep);
    save tmp_res
end

%plot(wts([1 length(wts)]),[opt_lik opt_lik],wts([1 ...
%      length(wts)]),[base_lik base_lik],wts(2:length(wts)),sdp_liks(2:end),'k');
%plot(wts([1 length(wts)]),[opt_lik opt_lik],wts,sdp_liks(1:end),'k.');
lw = 5;
fs = 14;
loglog(wts([1 length(wts)]),[opt_lik opt_lik],wts(1:end),sdp_liks(1:end),'k','LineWidth',lw);
set(gca,'FontSize',fs);
xlabel('\lambda','FontSize',fs);
ylabel('Likelihood','FontSize',fs);
print -depsc c:\proj\logemb\tex\Figures\trace_res_small.eps