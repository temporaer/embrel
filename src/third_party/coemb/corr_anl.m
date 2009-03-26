function [F,G,Du] = corr_anl(nxy,dim)

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

Q=P./sqrt(r*c);clear P

% [U,Dalpha,V] = svds(Q,dim);
[U,Dalpha,V] = svds(Q,dim);

%% Discard first elemt
A = diag(sqrt(r))*U;
B = diag(sqrt(c))*V;
Du = Dalpha;

for(i=1:nx), A1(i,:) =   A(i,:)/r(i);end
F = A1*Du(1:dim,:);
for(i=1:ny), B1(i,:) =   B(i,:)/c(i);end
G = B1*Du(1:dim,:);


F = F(:,2:end);
G = G(:,2:end);
dd = min([nx ny dim]);
Du = Du(2:dd,2:dd);

return

