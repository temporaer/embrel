function [new_phi,new_psi ] = red_const_from_edm(phi,psi,add_del)

NX = size(phi,1);
NY = size(psi,1);

k = NX+NY;
V = eye(k)-1/k*ones(k,1)*ones(1,k);

gram = [phi;psi]*[phi;psi]';

orig_dist = g_to_d(gram);
curr_dist = orig_dist;

add_const = 0;

%addmat = ones(size(curr_dist));
%addmat = addmat-eye(length(addmat));
addmat = zeros(size(curr_dist));
addmat(1:NX,NX+1:end) = 1;
addmat(NX+1:end,1:NX) = 1;
%addmat = ones(size(currdist));
%addmat(NX+1:end,NX+1:end) = 0;
%addmat(1:NX+1:end,NX+1:end) = 0;

cent_gram = -V*curr_dist*V;
%v = eigs_r11(gram,eye(size(cent_gram)),1);
%thr = 1e-1*v;

add_const = min(orig_dist(1:NX,NX+1:end));
add_const = 0.001*min(add_const);
for i=1
%    curr_dist(1:NX,NX+1:end) = orig_dist(1:NX,NX+1:end) - add_const;
%    curr_dist(NX+1:end,1:NX) = orig_dist(NX+1:end,1:NX) - add_const;
    curr_dist = orig_dist - addmat*add_const;
    cent_gram = -V*curr_dist*V;
    v = eigs_r11(cent_gram,eye(size(cent_gram)),7)
%    if min(real(v))<-thr
%        add_const = add_const-add_del;
%        break;
%    end
    add_const = add_const + add_del
end

%curr_dist(1:NX,NX+1:end) = orig_dist(1:NX,NX+1:end) - add_const;
%curr_dist(NX+1:end,1:NX) = orig_dist(NX+1:end,1:NX) - add_const;
%curr_dist = orig_dist - addmat*add_const;

%cent_gram = -V*curr_dist*V;

curr_dim = 2;
[v,d]=eig(cent_gram); %s_r11(cent_gram,eye(size(cent_gram)),15);
[v,d] = sortem(v,d);
% Generate reduced Gram matrix
newpoints = v(:,1:curr_dim)*sqrt(d(1:curr_dim,1:curr_dim));
new_phi = newpoints(1:NX,:);
new_psi = newpoints(NX+1:end,:);
