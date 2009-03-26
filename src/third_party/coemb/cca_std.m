function [F,G,rx,ry]=cca_std(X,Y)

% X: (n,dim_x)
% Y: (n,dim_y)
n = size(X,2);
sx = size(X,1);
sy = size(Y,1);

% Remove means
X = X - repmat(mean(X,2),1,n);
Y = Y - repmat(mean(Y,2),1,n);


SX = X*X'/n;% +  eps*eye(sx);
SY = Y*Y'/n;% + eps*eye(sy);

SXY = X*Y'/n;

A = zeros(sx+sy);
A(1:sx,sx+1:end) = SXY;
A(sx+1:end,1:sx) = SXY';

B = blkdiag(SX,SY);

[W,V] = eig(A,B);

WX = W(1:sx,1:sx);
WY = W(sx+1:end,sx+1:end);
F = WX'*X;
G = WY'*Y;
