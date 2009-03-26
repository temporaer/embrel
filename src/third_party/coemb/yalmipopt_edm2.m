function [final_gram,final_d,final_z] = yalmipopt_edm2(alld,NX,NY)
% A = mat_completion(Ainit,indices,C);
% dxy - Given distances between x and y. Must be kept up to a constant. All others are free

% There are NX^2+NY^2+1 variables, for the Gram of X,Y and the shift variable

dxy = alld(1:NX,NX+1:end);
dxx = alld(1:NX,1:NX);
dyy = alld(NX+1:end,NX+1:end);

N = NX+NY;

x0 = full([get_upper(dxx) get_upper(dyy)]');                         
ops = sdpsettings('solver','dsdp');
ops.dsdp.maxit = 5000;
%ops.csdp.maxiter = 1000;
%ops.csdp.minstepfrac = 0.7;
%ops.csdp.maxstepfrac = 1;
ops.csdp.usexzgap = 1;

P = sdpvar(N,N);
Z = sdpvar(1);
% P is the full Gram matrix. The constraints on it are:
% a. It should obey the relevant distances in the original matrix plus a constant
% b. It should be centered
% c. It should be PSD


F = set(P>9);      % Centering
%F = F + set(sum(P(:)==0));
%F = F +set(P(1,1)==0);
%F = F + set(Z==0);

if 1
for xi=1:NX
    for yi=1:NY
        fprintf('xi=%d yi=%d\n',xi,yi);
        F = F + set((P(xi,xi)+P(NX+yi,NX+yi)-P(xi,NX+yi)-P(NX+yi,xi))==(dxy(xi,yi)+Z));
    end
end
end

if 0
for xi1=1:NX+NY
    for xi2=xi1+1:NX+NY
        fprintf('xi1=%d xi2=%d\n',xi1,xi2);
        F = F + set((P(xi1,xi1)+P(xi2,xi2)-P(xi1,xi2)-P(xi2,xi1))==(alld(xi1,xi2)));
    end
end
end

sol = solvesdp(F,Z,ops);

final_gram = double(P);
final_z = double(Z);
final_d = g_to_d(final_gram);