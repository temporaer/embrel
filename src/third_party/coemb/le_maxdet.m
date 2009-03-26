function G = mat_completion(nxy)
% Find exponent of -D where D is the Euclidean distance Matrix

NX = size(nxy,1);
NY = size(nxy,2);

% check input arguments
n = NX+NY;

% Genereate matrices defining F, the exponent of the square Euclidean distances
F0 = eye(n); % Diagonal elements fixed to 1
F(:,1) = F0(:);

% Specify which elements are changed in the matrix
for i=1:NX+NY
    for j=i+1:NX+NY
        mat = zeros(NX+NY,NX+NY);
        mat(i,j) = 1;
        mat(j,i) = 1;
        F(:,end+1) = mat(:);
    end
end
F_blkszs = n;

nvar = size(F,2)-1;  % Decrease one for the constant factor

% Generate G. The diagonal matrix whose determinant we wish to maximize
% G(1,1) = sum(F(1:NX,NX+1:end))
NSamp = sum(nxy(:))
ng = 2*NSamp;

G0 = zeros(1,ng);
G(:,1) = G0(:);
dgind = NSamp+1;

for i=1:NX+NY
    for j=i+1:NX+NY 
        % Only need to specify the diagonal
        dg = zeros(1,ng);
        if i<=NX & j>NX
%            dg = zeros(1,ng);
            dg(1:NSamp) = 1;
%            mat = mat + diag(dg);
%            dg = zeros(1,ng);            
            dg(dgind:dgind+nxy(i,j-NX)-1)=1;
%            mat = mat + diag(dg);
            dgind = dgind + nxy(i,j-NX);
        end
        G(:,end+1) = dg'; % mat(:);
    end
end
G_blkszs = ones(1,ng);

% form c (constant part of the objective is ignored in the optimization)
c = zeros(nvar,1);                       % c_i = Tr G_i*C

gam = 100;
abstol = 1e-8;
reltol=0.5;
NIters = 1000;
[x,z,W,ul,infostr] = phase1(F,F_blkszs,G,G_blkszs,gam,abstol,reltol,NIters);
disp(infostr);

% call maxdet
x0 = x;                            % x0 = 0 is feasible
[x,Z,W,ul,hist,infostr] = maxdet(F,F_blkszs,G,G_blkszs,c,x0,...
      zeros(size(F,1),1),zeros(size(G,1),1),1e-8,0.5,100,100);

disp(infostr)

% form the solution D
D = reshape(F*[1;x],n,n);
%D = -log(D);
% Likelihood
pxy = D(1:NX,NX+1:end);
pxy = pxy/sum(pxy(:));

lik = sum(sum(nxy.*log(pxy)));
lik = lik/NSamp


