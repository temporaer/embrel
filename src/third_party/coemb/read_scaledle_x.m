function [phi,psi,beta] = read_scaledle_x(x,cg_params)

dim = cg_params.dim;
NX = cg_params.NX;
NY = cg_params.NY;

phi = reshape(x(1:NX*dim),NX,dim);
psi = reshape(x(NX*dim+1:end-1),NY,dim);
beta = x(end);
