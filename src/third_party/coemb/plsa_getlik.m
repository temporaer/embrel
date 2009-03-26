function lik = plsa_getlik(pwgz,pdgz,pz,data)

NZ = length(pz);
[d_tr,w_tr,v_tr] = find(data);
NTRAIN = length(v_tr);

pd = pz'*pdgz;

pzgdw = repmat(pz,1,NTRAIN).*pdgz(:,d_tr).*pwgz(:,w_tr);
lik = sum(pzgdw);

lik = lik./pd(d_tr);
