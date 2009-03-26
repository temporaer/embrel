function c =  make_c(pxy0,NX,NY,N)


n=N+1;
x = -1/(n+sqrt(n));
y = -1/sqrt(n);
a_diff = 2*(x-y)^2;
b_diff = 2*(x+1-y)*(x-y);
a_eq = (x-y)^2;
b_eq = (x+1-y)^2;

c = zeros(N,N);
dy = pxy0;

c(1:NX-1,NX:end) = -2*dy(2:end,:);
add = ones(NX-1,NY)*sum(dy(1,:))*a_diff+repmat(dy(1,:),NX-1,1)*2*(x-y);
c(1:NX-1,NX:end) = c(1:NX-1,NX:end) + add;

c(1:NX-1,1:NX-1) = sum(dy(1,:))*a_diff;

constx = (a_eq-a_diff)*sum(dy(1,:));
consty = a_eq*sum(dy(1,:));

diagX = sum(dy(2:end,:),2)+constx;

c = c + diag([diagX' zeros(1,NY)]);

c(NX:end,NX:end) = sum(dy(1,:))*a_diff;
diagY = sum(dy(2:end,:),1)+constx;
diagY = diagY+(b_eq - a_eq)*dy(1,:);
c = c + diag([zeros(1,NX-1) diagY]);

v = dy(1,1:NY);
e = ones(NY-1,1);
k = (v(1:end-1)'*e' + e*v(2:end))*2*(x-y);

c(NX:end-1,NX+1:end) = c(NX:end-1,NX+1:end) + triu(k);
