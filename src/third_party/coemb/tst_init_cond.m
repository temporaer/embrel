% Generate joint
rand('seed',0);
NX = 20;
NY = 35;
nxy = rand(NX,NY);
nxy = nxy/sum(nxy(:));

params.noise = 1e-12;
params.x_marg = 'F';
params.y_marg = 'M';
params.b_cond = 1;
pxy0s{1} = nxy;
dim_list = 2;

[phi,psi,as,bs,lik] = code(pxy0s,dim_list,params);

d = my_dist(phi{2},psi{2}');