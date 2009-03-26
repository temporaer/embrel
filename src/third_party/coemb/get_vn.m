function vn = get_vn(n)

%vn = 1/sqrt(2)*[-1*ones(1,n-1);eye(n-1)];
% Following Alfakih
x = -1/(n+sqrt(n));
y = -1/sqrt(n);
xm = speye(n-1);
xm = xm+x;
vn =[(y*ones(1,n-1));xm];
