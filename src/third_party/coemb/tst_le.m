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
xv = rand(NX,dimension);
yv = rand(NY,dimension);
%sc=1e-1;
%xv = xv+rand(size(xv))*sc;
%yv = yv+rand(size(yv))*sc;

real_d = dist(xv,yv'); 
% Generate distribution
pxy = exp(-real_d.^2);
pxy = pxy/sum(pxy(:));
% Generate data

linsamp = draw_discrete_c(pxy,N);
[x_samp,y_samp]  = ind2sub([NX,NY],linsamp);
samp = [x_samp;y_samp];

t = lecent_solve_beta2(samp,NX,NY);

Z = sum(pxy(:));
logpxy = log(pxy); %-curr_dmat-log(Z);  %pxy/sum(pxy(:));
loglik = 0;
for i=1:length(samp)
   loglik = loglik + logpxy(samp(1,i),samp(2,i));
end
loglik = loglik/length(samp)

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


