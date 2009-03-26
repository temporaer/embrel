function [phi,psi,dist_mat] = run_le_proj(nxy,curr_dim,phi0,psi0)

% Run in full dimension and project

NX = size(nxy,1);
NY = size(nxy,2);
pxy0 = nxy/sum(nxy(:));
N = NX+NY;

if exist('phi0')
    xv0 = phi0;
    yv0 = psi0;
else
    xv0 = rand(NX,N)*1e-3;
    yv0 = rand(NY,N)*1e-3;
end

x0 = [xv0(:);yv0(:)];

cg_params.NX = NX;
cg_params.NY = NY;
cg_params.dim = N;
cg_params.pxy0 = pxy0;

opts = -1*zeros(1,9);
x = conj_grad('logemb_grad',cg_params,x0,opts);
[phi,psi] = read_legrad_x(x,cg_params);

%phi = phi-repmat(mean(phi),NX,1);
%psi = psi-repmat(mean(psi),NY,1);
% Generate Gram matrix
t = [phi;psi]*[phi;psi]';

dist_mat = g_to_d(t);

[v,d] = eig(t);
[v,d] = sortem(v,d);
diag(d)

% Generate reduced Gram matrix
newpoints = v(:,1:curr_dim)*sqrt(d(1:curr_dim,1:curr_dim));
phi = newpoints(1:NX,:);
psi = newpoints(NX+1:end,:);


