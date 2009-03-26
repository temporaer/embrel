function [cx,cxeq,g_const,g_const_eq] = marg_model_constr(x,params)

g_const = [];
dim = params.dim;
NX = params.NX;
NY = params.NY;
pxy0 = params.pxy0;
px0 = params.px0;
py0 = params.py0;
oprod = params.oprod;

%[v,j] = max(pxy0(:));
%[i,j] = ind2sub(size(pxy0),j);

[phi,psi] = read_legrad_x(x,params);
phi_tr = phi';
psi_tr = psi';

% Generate model distribution
%d = my_dist(phi,psi_tr);
%mx = max(d(:))*mybeta;
%logpxy = mx-mybeta*d-log(sum(sum(oprod.*exp(-d*mybeta+mx))))+log(oprod); %+log(px0*py0);
[f,gphi,gpsi,ga,gb,logpxy] = mrgle_grad_fast_prm(phi,psi,log(px0),log(py0),pxy0,px0,py0,0);
pxy = exp(logpxy);

px = sum(pxy,2);
py = sum(pxy,1);

cx = []; %zeros(NX+NY,1)-1;        % Inequality constraint. Set to minus one so that it is always true
cxeq = [px-px0;(py-py0)'];

if 1
psi_exp_d = make_cond_dist(pxy,0)*psi;
phi_exp_d = (phi_tr*make_cond_dist(pxy,1))';

g_const_eq=[];
g_const_eq = zeros(NX+NY,(NX+NY)*dim);
ind = 1;

for xi=1:NX
    % Derivative of p(x) w.r.t \phi(x)
    v = ones(NX,1)*px(xi);
    v(xi) = v(xi) - 1;
    dpx_phix = (phi-psi_exp_d).*repmat(v,1,dim).*repmat(px,1,dim);
    
    % Derivative of p(x) w.r.t \psi(y)
    
    dpx_psiy = px(xi)*(psi-phi_exp_d).*repmat(py',1,dim)-repmat(pxy(xi,:)',1,dim).*(psi-repmat(phi(xi,:),NY,1));
    g_const_eq(ind,:) = [dpx_phix(:); dpx_psiy(:)]';
    ind = ind+1;
end

for yi=1:NY
    % Derivative of p(y) w.r.t \phi(x)    
    v = ones(NY,1)*py(yi);
    v(yi) = v(yi) - 1;
    dpy_psiy = (psi-phi_exp_d).*repmat(v,1,dim).*repmat(py',1,dim);
    
    % Derivative of p(y) w.r.t \psi(x)    
    dpy_phix = py(yi)*(phi-psi_exp_d).*repmat(px,1,dim)-repmat(pxy(:,yi),1,dim).*(phi-repmat(psi(yi,:),NX,1));
    
    g_const_eq(ind,:) = [dpy_phix(:); dpy_psiy(:)]';    
    ind = ind + 1;
end

g_const_eq = 2*g_const_eq';

end


return;
px0 = sum(pxy0,2);
py0 = sum(pxy0,1);

gphi = diag(px)*phi-pxy*psi - (diag(px0)*phi-pxy0*psi);
%gpsi2 = diag(py)*psi-pxy'*phi - (diag(py0)*psi-pxy0'*phi);
gpsi = psi_tr*diag(py)-phi_tr*pxy - (psi_tr*diag(py0)-phi_tr*pxy0);
gpsi = gpsi';

g = [-gphi(:);-gpsi(:)];

%f = log(Z)+sum(pxy0(:).*d(:));
f = -sum(pxy0(:).*logpxy(:))+params.lik0;




