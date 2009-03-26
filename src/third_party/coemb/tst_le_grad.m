%seed
%rand('seed',seed)
% Generate data - Embed X and Y into  a circle in R2
N = 2000;
NX = 10;
NY = 10;
% Generaete embeddings
for i=1:NX
    [xv(i,1) xv(i,2)] = pol2cart((2*pi)/NX*(i-1),1);
end
for i=1:NY
    [yv(i,1) yv(i,2)]  = pol2cart((2*pi)/NY*(i-1),1);
end
dimension = 2;
%xv = rand(NX,dimension);
%yv = rand(NY,dimension);
%sc=1e-1;
%xv = xv+rand(size(xv))*sc;
%yv = yv+rand(size(yv))*sc;

real_d = dist(xv,yv'); 
% Generate distribution
pxy = exp(-real_d.^2);
pxy = pxy/sum(pxy(:));
% Generate data

linsamp = draw_from_p(pxy(:),N);
[x_samp,y_samp]  = ind2sub([NX,NY],linsamp);
samp = [x_samp;y_samp];
% Generate empirical sample
for i=1:NX
    for j=1:NY
        pxy0(i,j) = sum(samp(1,:)==i & samp(2,:)==j);
    end
end
pxy0 = pxy0/sum(pxy0(:));
pxy0 = pxy


max_repeats=5;
for curr_dim = 2
  for(i_repeat=1:max_repeats)
    fprintf('i_repeat=%d\n',i_repeat);

    xv0 = rand(NX,curr_dim);
    yv0 = rand(NY,curr_dim);
  
    x0 = [xv0(:);yv0(:)];
    cg_params.NX = NX;
    cg_params.NY = NY;
    cg_params.dim = curr_dim;
    cg_params.pxy0 = pxy0;
    
    opts = -1*zeros(1,9);
    x = conj_grad('logemb_grad',cg_params,x0,opts);
    [phis{i_repeat},psis{i_repeat}] = read_legrad_x(x,cg_params);
    liks(i_repeat) = -logemb_grad(x,cg_params);
  end
  [max_lik,bst_repeat] = max(liks);
  lik(curr_dim) = max_lik;
  phi = phis{bst_repeat};
  psi = psis{bst_repeat};  
end
% liks(seed) = lik(2);

clf
subplot(2,1,1);
plot(phi(:,1),phi(:,2),'.');
hold on
plot(psi(:,1),psi(:,2),'r.');
axis equal
subplot(2,1,2);
clear i;
plot(angle(phi(:,1)+i*phi(:,2)),angle(psi(:,1)+i*psi(:,2)),'.');
axis tight
return;
Z = sum(pxy(:));
logpxy = log(pxy); %-curr_dmat-log(Z);  %pxy/sum(pxy(:));
loglik = 0;
for i=1:length(samp)
   loglik = loglik + logpxy(samp(1,i),samp(2,i));
end
loglik = loglik/length(samp)

t = [phi;psi]*[phi;psi]';
d = g_to_d(t);
subd = d(1:NX,NX+1:end);
found_pxy = exp(-subd);
%found_pxy = found_pxy(1:NX,NX+1:end);
found_pxy = found_pxy/sum(found_pxy(:));
logpxy = log(found_pxy); %-curr_dmat-log(Z);  %pxy/sum(pxy(:));
found_loglik = 0;
for i=1:length(samp)
   found_loglik = found_loglik + logpxy(samp(1,i),samp(2,i));
end
found_loglik = found_loglik/length(samp)

