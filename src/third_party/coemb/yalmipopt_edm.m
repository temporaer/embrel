function final_gram = opt_edm(alld,NX,NY)
% A = mat_completion(Ainit,indices,C);
% dxy - Given distances between x and y. Must be kept up to a constant. All others are free

% There are NX^2+NY^2+1 variables, for the Gram of X,Y and the shift variable

dxy = alld(1:NX,NX+1:end);
dxx = alld(1:NX,1:NX);
dyy = alld(NX+1:end,NX+1:end);

N = NX+NY;


G=[];

for x1=1:NX
    for x2=x1+1:NX
        m = zeros(N,N);
        m(x1,x2)=1;
        m(x2,x1)=1;
        m = m  + (2/N^2)*ones(N,N);
        tmp = zeros(N,N);
        tmp(x1,:) = 1/N;
        tmp(:,x1) = tmp(:,x1)+1/N;
        m = m - tmp;
        tmp = zeros(N,N);
        tmp(:,x2) = 1/N;
        tmp(x2,:) = tmp(x2,:)+1/N;
        m = m - tmp;
        G(end+1,:) = m(:)';
    end
end
    
for y1=1:NY
    for y2=y1+1:NY
        m = zeros(N,N);
        m(NX+y1,NX+y2)=1;
        m(NX+y2,NX+y1)=1;
        m = m + (2/N^2)*ones(N,N);
        tmp = zeros(N,N);
        tmp(NX+y1,:) = 1/N;
        tmp(:,NX+y1) = tmp(:,NX+y1) + 1/N;
        m = m - tmp;
        tmp = zeros(N,N);
        tmp(:,NX+y2) = 1/N;
        tmp(NX+y2,:) = tmp(NX+y2,:)+1/N;
        m = m - tmp;
        G(end+1,:) = m(:)';
    end
end

G0 = 0;
G2 = 0;
% Collect all the free elements
for xi=1:NX
    for yi=NX+1:N
        m = zeros(N,N);
        m(xi,yi) = 1;
        m(yi,xi) = 1;
        m = m +(2/N^2)*ones(N,N);
        tmp = zeros(N,N);
        tmp(xi,:) = 1/N;
        tmp(:,xi) = tmp(:,xi) + 1/N;
        m = m - tmp;
        tmp = zeros(N,N);
        tmp(:,yi) = 1/N;
        tmp(yi,:) = tmp(yi,:) + 1/N;
        m = m - tmp;
        G0  = G0 + dxy(xi,yi-NX)*m(:)';
        G2 = G2+m(:)';
    end
end

% Add G0 and G for Z variable
G = [G0;G;G2];
%G = [G0;G;];
G = -G';     
G_blkszs = N;
NVars = size(G,2)-1;

x0 = full([get_upper(dxx) get_upper(dyy)]');                         

x = sdpvar(NVars,1);
P = reshape(G*[1;x],N,N);

F = set(P>0);
F = F + set(sum(P(:))==0);
sol = solvesdp(F,x(end));

final_gram = reshape(G*[1;double(x)],N,N);
