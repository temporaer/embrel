function [f,g,basef,newf] = sdp_grad(x,params)

%A = params.A;
%c = params.c;
trace_grad = params.trace_grad;
%F = params.F;
nxy = params.nxy;
tr_wgt = params.tr_wgt;
t = params.t;
uind = params.uind;
NX = params.NX;
NY = params.NY;
c = params.c;
N = NX+NY;

%y = A*x;
%gradf = exp(y-y(1))./sum(exp(y-y(1)));
% Create distance matrix in newmat
gram_mat = vec_to_smat_fast(x,N); %F*x;  %
%N = sqrt(size(F,1));
gram_mat = reshape(gram_mat,N,N);

ev = eig(gram_mat);

if (min(ev)<=0)
    f = Inf;
    g = zeros(size(x));
    return;
end

gram_diag = diag(gram_mat);
dist_mat = gram_diag*ones(1,N) + ones(N,1)*gram_diag' - 2*gram_mat;
% Dist mat reduced to [NX,NY]
%red_dist_mat = dist_mat(1:NX,NX+1:end);
red_dist_mat2 = dist_mat(NX+1:end,1:NX);
% Create the gram matrix
  %    basef = log(sum(exp(y-y(1))))+y(1)-dot(c,x); 
%    logZ = log(sum(exp(-red_dist_mat(:)+red_dist_mat(1,1))))-red_dist_mat(1,1);
logZ = log(sum(exp(-red_dist_mat2(:)+red_dist_mat2(1,1))))-red_dist_mat2(1,1);    
basef = logZ+dot(x,c);
%    sum(red_dist_mat(:).*nxy(:));  Should be like dot(x,c)
newf =  basef + tr_wgt*trace(gram_mat);    
f = t*newf - sum(log(ev));    %det(gram_mat));     

if nargout==1
  return;
end


bar_grad = det_barrier_grad_gd(gram_mat);

%g = A' * (gradf*t)-t*c +bar_grad + trace_grad; 

% Convert to gradient 
if 0 
Zgrad = sparse(size(dist_mat,1),size(dist_mat,2));
Zgrad(1:NX,NX+1:end) = exp(-red_dist_mat);
xsum = sum(Zgrad,2);
ysum = sum(Zgrad,1);
Zgrad = Zgrad - diag(xsum);
Zgrad = Zgrad - diag(ysum);
Zgrad(1:NX,NX+1:end) = 2*Zgrad(1:NX,NX+1:end);
Zgrad(NX+1:end,1:NX) = Zgrad(1:NX,NX+1:end)';
%Zgrad = Zgrad(uind)/exp(logZ);
Zgrad = smat_to_vec_c(full(Zgrad))/exp(logZ);
g = t*(Zgrad'+c'+trace_grad)+ bar_grad  ; 
end

Zgrad2 = sparse(size(dist_mat,1),size(dist_mat,2));
Zgrad2(NX+1:end,1:NX) = exp(-red_dist_mat2-logZ);
xsum = sum(Zgrad2,2);
ysum = sum(Zgrad2,1);
Zgrad2 = Zgrad2 - diag(xsum);
Zgrad2 = Zgrad2 - diag(ysum);
Zgrad2(NX+1:end,1:NX) = 2*Zgrad2(NX+1:end,1:NX);
Zgrad2 = smat_to_vec_c(full(Zgrad2));

g = t*(Zgrad2'+c'+trace_grad)+ bar_grad  ; 

%if isnan(f) | sum(isnan(g))>1
%   disp('Funny numbers');
%   keyboard;     
%end
