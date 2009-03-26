% Transform Gram matrix X*X' to distance matrix d(i,j)=|x_i-x_j|^2

function d = g_to_d(g)

N = size(g,1);
dg = diag(g);

d = repmat(dg',N,1) + repmat(dg,1,N) - 2*g;

% d = zeros(size(g));
% for i=1:size(d,1)
%     for j=1:size(d,2)
%         d(i,j) = g(i,i)+g(j,j)-2*g(i,j);
%     end
% end