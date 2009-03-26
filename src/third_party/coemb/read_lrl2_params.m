function [phi,psi,c] = read_lrl2_params(x,cg_params)

dim = cg_params.dim;
x_i = 1;

% Read the phi, which is common to all distributions
sz = cg_params.NX*dim;
phi = reshape(x(x_i:x_i+sz-1),cg_params.NX,dim);
x_i = x_i+sz;
sz = cg_params.NY*dim;
psi = reshape(x(x_i:x_i+sz-1),cg_params.NY,dim);
c = x(end);
