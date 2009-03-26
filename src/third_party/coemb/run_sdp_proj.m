function [phi,psi,ev] = run_sdp_proj(nxy,curr_dim,tr_wgt)

% Run in full dimension and project

NX = size(nxy,1);
NY = size(nxy,2);
pxy0 = nxy/sum(nxy(:));
N = NX+NY;

% xv0 = rand(NX,N)*1e-3;
% yv0 = rand(NY,N)*1e-3;
% 
% x0 = [xv0(:);yv0(:)];
% Returns the Gram matrix
%g = ler_solve(nxy,NX,NY);
%[phi0,psi0] = logemb_grad_dim(pxy0,curr_dim,1);
%g0 = [phi0;psi0]*[phi0;psi0]';
%e0 = eig(g0);
%if min(e0)<=0
%    g0 = g0 - 2*eye(length(g0))*(min(e0)-eps);
%end

g = le_trace_gd_str(nxy,NX,NY,tr_wgt);
%%g = le_trace_solve_gd(nxy,NX,NY,tr_wgt);
%train_lik = 0;
%[g] = le_hess(nxy,NX,NY,tr_wgt);
%g = le_trace_solve_fixorig(nxy,NX,NY,tr_wgt);

%g = lecent_solve(nxy,NX,NY);

dist_mat = g_to_d(g);
% Center the gram matrix
%g = g-sum(g(:));

[v,d] = eig(g);
[v,d] = sortem(v,d);
ev = diag(d);

% Generate reduced Gram matrix
newpoints = v(:,1:curr_dim)*sqrt(d(1:curr_dim,1:curr_dim));
phi = newpoints(1:NX,:);
psi = newpoints(NX+1:end,:);


