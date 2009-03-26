function [phi,psi] = reduce_gram(gram,curr_dim,NX,NY)


[v,d] = eig(gram);
[v,d] = sortem(v,d);


% Generate reduced Gram matrix
newpoints = v(:,1:curr_dim)*sqrt(d(1:curr_dim,1:curr_dim));
phi = newpoints(1:NX,:);
psi = newpoints(NX+1:end,:);

