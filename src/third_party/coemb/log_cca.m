function [F,G,Du] = log_cca(nxy,dim,option)

% Perform correspondance analysis on give count matrix
nx = size(nxy,1);
ny = size(nxy,2);

% Joint Distribution
thr = 1e-12;
nxy = nxy+thr;
P = nxy/sum(nxy(:));clear nxy

r = sum(P,2); % px
c = sum(P,1); % py


if(~isempty(find(c==0)))
  error('corr_anl: nxy contains zero columns.');
end

P = log(P);
lambda = P-repmat(mean(P,2),1,ny)-repmat(mean(P,1),nx,1)+mean(P(:));

% [U,Dalpha,V] = svds(Q,dim);
% Notation as in Cox and Cox
[A,Dalpha,B] = svds(lambda,dim);

% Remove means

% Make A orthogonal w.r.t. inv(diag(r))
%A = diag(1./sqrt(r))*A;
%B = diag(1./sqrt(c))*B;

F = A;
G = B;
Du = 0;

return

