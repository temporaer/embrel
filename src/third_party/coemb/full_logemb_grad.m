function [f,g] = full_logemb_grad(x,params) 

dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;
pxx0 = params.pxx0;
pyy0 = params.pyy0;
[phi,psi] = read_legrad_x(x,params);
% Get pxx0 gradient

if ~isempty(pxx0)
	new_params.dim = dim;
	new_params.NX = NX;
	new_params.NY = NX;
	new_params.pxy0 = pxx0;
	new_x = [phi(:);phi(:)];
	[pxx_val,pxx_grad] = logemb_grad(new_x(:),new_params);
else
    pxx_val = 0;
    pxx_grad = zeros(1,length(phi(:))*2);
end

if ~isempty(pyy0)
	new_params.NX = NY;
	new_params.NY = NY;
	new_params.pxy0 = pyy0;
	new_x = [psi(:);psi(:)];
	[pyy_val,pyy_grad] = logemb_grad(new_x(:),new_params);
else
    pyy_val = 0;
    pyy_grad = zeros(1,length(psi(:))*2);
end

new_params.NX = NX;
new_params.NY = NY;
new_params.pxy0 = pxy0;
[pxy_val,pxy_grad] = logemb_grad(x,new_params);

f = params.cxx*pxx_val + params.cyy*pyy_val + params.cxy*pxy_val;
g = zeros(size(x));
g(1:NX*dim) = params.cxx*pxx_grad(1:NX*dim);
g(NX*dim+1:end) = params.cyy*pyy_grad(1:NY*dim);
g = g + params.cxy*pxy_grad;
%g = params.cxx*pxx_grad + params.cyy*pyy_grad + params.cxy*pxy_grad;
