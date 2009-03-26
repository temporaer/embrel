function [F,G,Du] = ca_cca(nxy,dim,option)

% Due to a constant factor in the SVD
dim = dim+1;

% Perform correspondance analysis on give count matrix
nx = size(nxy,1);
ny = size(nxy,2);

% Joint Distribution
P = nxy/sum(nxy(:));clear nxy

r = sum(P,2); % px
if(~isempty(find(r==0)))
  error('corr_anl: nxy contains zero columns.');
end
c = sum(P,1); % py
if(~isempty(find(c==0)))
  error('corr_anl: nxy contains zero columns.');
end

Q=P./sqrt(r*c);

% [U,Dalpha,V] = svds(Q,dim);
% Notation as in Cox and Cox
[A,Dalpha,B] = svds(Q,dim);

% Make A orthogonal w.r.t. inv(diag(r))
A = diag(1./sqrt(r))*A;
B = diag(1./sqrt(c))*B;

F = A(:,2:end);
G = B(:,2:end);
dd = min([nx ny dim]);
Du  = Dalpha;
return

