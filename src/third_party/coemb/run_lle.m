function [F,G,lik] = run_lle(nxy,dim,params)

% [Y] = lle(X,K,dmax)
%
% X = data as D x N matrix (D = dimensionality, N = #points)
% K = number of neighbors
% dmax = max embedding dimensionality
% Y = embedding as dmax x N matrix

% Perform correspondance analysis on give count matrix
nx = size(nxy,1);
ny = size(nxy,2);

nxy = nxy+rand(size(nxy))*params.tol;
% Return an embedding of Y, which is the documents in our case
G  = lle(nxy,params.k,dim,params.b_l2dist);
G = G';

F = zeros(nx,dim);
lik = 0;
