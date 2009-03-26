function [f,g]=lrl2_grad(x,params)

[phi,psi,c] = read_lrl2_params(x,params);
NX = size(phi,1);
NY = size(psi,1);
N = NX*NY;

lrat = params.lrat;

d = my_dist(phi,psi');

% Here we are modeling -log p(x,y)/(p(x)p(y)) as d(\phi,\psi)-c
if params.b_log
  delta = lrat+c-d;
  
  fx = sum(delta,2);
  fy = sum(delta,1);
  
  mfx = spdiags(fx,0,NX,NX);
  mfy = spdiags(fy',0,NY,NY);
  
  grad_phi = 2*(-mfx*phi+delta*psi);
  grad_psi = 2*(-psi'*mfy+phi'*delta);
  grad_psi = grad_psi';
  grad_c = sum(delta(:));
else
  % Here we are modeling  p(x,y)/(p(x)p(y)) as exp(-d(\phi,\psi)+c)
  expd = exp(c-d);
  delta = lrat-expd;
  delta2 = delta.*expd;
  
  fx = sum(delta2,2);
  fy = sum(delta2,1);

  mfx = spdiags(fx,0,NX,NX);
  mfy = spdiags(fy',0,NY,NY);
  
  grad_phi = -2*(-mfx*phi+delta2*psi);
  grad_psi = -2*(-psi'*mfy+phi'*delta2);
  grad_psi = grad_psi';
  grad_c = -sum(delta2(:));
end

g = [grad_phi(:); grad_psi(:);grad_c]/N;

f = sum(delta(:).^2)/N;

