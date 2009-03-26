function [lik_ds,w,phis,psis] = read_sdp_tmp(fn,ds)  %le_solve(samp,NX,NY)

load(fn,'x','py0','pxy0','big_iters','NX','NY','N','tr_wgt');
w = tr_wgt;
lik_ds = -10;


if big_iters < 4
  return;
end


% This is (n-1) SDP
m = vec_to_smat_fast(x(1:end-NX),N-1);
% Transform it into distance matrix
dist_m = kappav(m);
zxs = x(end-NX+1:end); 

clear m x


% Translate to centered Gram
k = NX+NY;
cent_m = dist_m -1/k*repmat(sum(dist_m,1),N,1);
cent_m = cent_m - 1/k*repmat(sum(cent_m,2),1,N);

%V = eye(k)-1/k*ones(k,1)*ones(1,k);

%cent_m = -0.5*V*dist_m;

%cent_m = cent_m*V;
%cent_m = -0.5*V*dist_m*V;

m = -cent_m;
clear cent_m;


%sq = sqrtm(full(cent_m));
%sq  = sq-repmat(mean(sq),size(sq,1),1);
%m = sq*sq'; %m-repmat(mean(m),size(m,1),1);


pxy = exp(-dist_m(1:NX,NX+1:end)-repmat(zxs,1,NY)+repmat(log(py0),NX,1));
xlse = full(sum(pxy,2)); 



print_to_log('wt=%g t=%d NormTo=%g-%g. Lost=%g\n',w,big_iters,min(xlse),max(xlse),mean(xlse<0.9));

opts.disp = 0;
[v,d] = eigs(m,max(ds),'LM',opts);
[v,d] = sortem(v,d);
evs = diag(d);

% Generate reduced Gram matrix
for di = 1:length(ds)
  curr_dim = ds(di);
  newpoints = v(:,1:curr_dim)*sqrt(d(1:curr_dim,1:curr_dim));
  phis{ds(di)} = real(newpoints(1:NX,:));
  psis{ds(di)} = real(newpoints(NX+1:end,:));
  b_ycond = 0;
  lik_ds(di) =full(get_le_cond_lik(phis{ds(di)},psis{ds(di)},pxy0,b_ycond));  
end

% Get its likelihood on the current data





