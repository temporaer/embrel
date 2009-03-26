function [phi,psi,a,b] = read_legrad_x(x,cg_params)

dim = cg_params.dim;
NX = cg_params.NX;
NY = cg_params.NY;

phi = reshape(x(1:NX*dim),NX,dim);
start = NX*dim+1;
psi = reshape(x(start:start+NY*dim-1),NY,dim);
start = start+NY*dim;
a   = x(start:start+NX-1);
start = start+NX;
b   = x(start:start+NY-1);
