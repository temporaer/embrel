function [F,G]=cca(X,Y,b_redmean)


% X: (n,dim_x)
% Y: (n,dim_y)
n = size(X,2);
sx = size(X,1);
sy = size(Y,1);

% Remove means
if nargin>2 & b_redmean
  X = X - repmat(mean(X,2),1,n);
  Y = Y - repmat(mean(Y,2),1,n);
end


SX = X*X'/n;% +  eps*eye(sx);
SY = Y*Y'/n;% + eps*eye(sy);

SXY = X*Y'/n;

cca_mat = inv(SX)*SXY*inv(SY)*SXY';
cca_mat = inv(sqrtm(SX))*SXY*inv(SY)*SXY'*inv(sqrtm(SX));

[Wx,rx,f] = svd(cca_mat);
[Wx,rx] = sortem(Wx,abs(rx));
F = Wx'*inv(sqrtm(SX));

cca_mat = inv(SY)*SXY'*inv(SX)*SXY;
cca_mat = inv(sqrtm(SY))*SXY'*inv(SX)*SXY*inv(sqrtm(SY));
[Wy,ry,f] = svd(cca_mat);
[Wy,ry] = sortem(Wy,abs(ry));
Wy = Wy'*inv(sqrtm(SY));

G = Wy;

return;












SX = X*X'/n;
SY = Y*Y'/n;

A = zeros(sx+sy);
A(1:sx,sx+1:end) = SXY;
A(sx+1:end,1:sx) = SXY';

B = blkdiag(SX,SY);

[F2,G2] = eig(A,B);
keyboard