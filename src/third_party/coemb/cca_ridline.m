function [F,G,rx,ry]=cca_ridline(X,Y)

% Get tid of the last rows of X and Y since they are known to be dependent
X = X(1:end-1,:);
Y = Y(1:end-1,:);

% Run standard cca on it
[F,G] = cca(X,Y,1);

